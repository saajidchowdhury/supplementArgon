clear all;
setenv('TZ', 'America/New_York');
fclose('all');
load("interpolation.mat");
main = string(datetime('now','Format','M-d-y@HH.mm.ss'))+"dmcread2B";
mkdir(main);
tstart = tic;


%dmc parameters
ns = [4:13]; %[4:58] number of particles in paper, spread across different runs
cores = 10; %28 or 60 cores in paper 
coresPerN = 1; %25 or more cores per cluster size in paper
% parpool('Processes',cores);                            %uncomment parpool
ln = length(ns);
runs = ln*coresPerN;
N0s = 200 * ones(ln,1); %20000 walkers in paper
dts = 10 * ones(ln,1);
t0s = 500 * dts; %t = 5000 a.u. of (imaginary) time in paper
maxlength = 1e8;


%basin hopping parameters
T = 0.001;
R = 100;
stepsize = 7.5;
niter = 100;
timefactor = N0s(1)/20000 * t0s(1)/dts(1)/15000 * 23.5 / 15.875;


%potential parameters
xfrozen = [0,0,0]; %in Bohr
Eh = 27.211386245981; %Hartree energy in eV
a0 = 0.529177210544; %Bohr radius in angstroms
me = 5.485799090441e-4; %electron mass in u
%atom-atom
epsilon = 1.23e-2 / Eh;
sigma = 3.357 / a0;
%short range
U = 4.3 / Eh;
re = 1.34 / a0;
a1 = 0.475; b1 = 1.370 * a0; %r < re
a2 = 0.400; b2 = 1.890 * a0; %r > re
%switching function on interval [3.6/a0, 3.86/a0] = [r1,r2]
r1 = 3.6 / a0;
r2 = 3.86 / a0;
A = -0.095 / Eh;
B = 0.01 / Eh * a0;
C = -0.085 / Eh;
gamma = 7 * a0;
rs = 3.5 / a0;
%long range VIII(r) = -alpha/2r^4
alpha = 11.08;
mp = 1836.152673426;
mAr = 39.9623831237 / me;


%job assignment
jobs = cell(runs,1);
runtimes = zeros(runs,1);
for i = 1:ln
    for j = 1:coresPerN
        runtimes((i-1)*coresPerN + j) = (39.9*ns(i)^2 + 500.9*ns(i) + 8177) * timefactor;
        jobs{(i-1)*coresPerN + j} = [ns(i),N0s(i),dts(i),t0s(i),i]; %attributes
    end
end
assignments = cell(cores,1);
totaltimes = zeros(cores,1);
[~,sortedIDs] = sort(runtimes,'descend');
for i = sortedIDs'
    [time,core] = min(totaltimes);
    assignments{core} = [jobs{i}; assignments{core}];
    totaltimes(core) = totaltimes(core) + runtimes(i);
end
for i = 1:cores
    fprintf("Core %d will do: ", i);
    for job = assignments{i}'
        fprintf("%d ",job(1));
    end
    fprintf(". Runtime = %.20f sec.\n", totaltimes(i));
end
maxtime = max(totaltimes);
fprintf("It will take %.20f seconds (%.20f hours).\n",maxtime,maxtime/3600);


%allocate and send jobs
data = zeros(ln,8); %attributes
dataPF = cell(cores,1);
parfor corei = 1:cores
    ass = assignments{corei};
    [jobn,~] = size(ass);
    for jobi = 1:jobn
        tp = tic;
        job = ass(jobi,:);
        n = job(1); N0 = job(2); dt = job(3); t0 = job(4); datai = job(end); %attributes
        sub = sprintf("n=%dcore=%djob=%d", n, corei, jobi);
        mkdir(main+"/"+sub);
        Nmax = 10*N0;
        steps = 2*floor(t0/dt);


        %reading
        fprintf("Core=%d n=%d N0=%d dt=%d t0=%d: Starting reading.\n",corei,n,N0,dt,t0);
        xsol = zeros(1,3*n);
        f = fopen(sprintf("nArNH+2B/Ar%dH+.xyz",n-1));
        nread = fscanf(f, "%d");
        assert(n == nread);
        fscanf(f, "%[H]");
        xsol(1:3) = fscanf(f, "%f");
        for i = 2:n
            fscanf(f, "%[Ar]");
            xsol(3*(i-1)+1:3*(i-1)+3) = fscanf(f, "%f");
        end
        Vsol = V(xsol,epsilon,sigma,U,re,a1,b1,a2,b2,r1,r2,A,B,C,gamma,rs,alpha);
        fprintf("Core=%d n=%d N0=%d dt=%d t0=%d: Finished reading. Vmin = %.20f.\n",corei,n,N0,dt,t0,Vsol);


        %dmc calculation
        ER = 0;
        ERlist = zeros(steps,1);
        ERlist(1) = ER;
        inters = zeros(maxlength,1);
        dists = zeros(maxlength,1);
        interi = 1;
        disti = 1;
        flag = zeros(Nmax, 1);
        x = zeros(Nmax, 3*n);
        for i = 1:N0
            flag(i) = 1;
            x(i,:) = xsol;
        end
        q = zeros(Nmax,1);
        for i = N0+1:Nmax
            q(i-N0) = i;
        end
        qi = Nmax-N0;
        for t = 1:steps
            for i = 1:Nmax
                if flag(i) ~= 0
                    flag(i) = 1;
                    x(i,1:3) = x(i,1:3) + sqrt(dt/mp)*randn(1,3);
                    x(i,4:3*n) = x(i,4:3*n) + sqrt(dt/mAr)*randn(1,3*n-4+1);
                end
            end
            N1 = 0;
            Vsum = 0;
            for i = 1:Nmax
                if flag(i) == 1
                    Vx = V(x(i,:),epsilon,sigma,U,re,a1,b1,a2,b2,r1,r2,A,B,C,gamma,rs,alpha);
                    W = exp(-(Vx-ER)*dt);
                    mn = min(floor(W+rand),3);
                    if mn == 0
                        flag(i) = 0;
                        q(qi+1) = i;
                        qi = qi+1;
                    elseif mn == 2
                        top = q(qi);
                        qi = qi-1;
                        flag(top) = 2;
                        x(top,:) = x(i,:);
                    elseif mn == 3
                        top1 = q(qi);
                        qi = qi-1;
                        top2 = q(qi);
                        qi = qi-1;
                        flag(top1) = 2;
                        x(top1,:) = x(i,:);
                        flag(top2) = 2;
                        x(top2,:) = x(i,:);
                    end
                    N1 = N1 + mn;
                    Vsum = Vsum + mn*Vx;
                end
            end
            assert(N1 ~= 0);
            V1 = Vsum/N1;
            ER = V1 + 1/dt*(1-N1/N0);
            ERlist(t) = ER;
            if t >= 0.75*steps
                for i = 1:Nmax
                    if flag(i) ~= 0
                        for j = 1:n-1
                            if disti <= maxlength
                                dists(disti) = sqrt((x(i,3*j+1)-x(i,3*0+1))^2 + ...
                                                    (x(i,3*j+2)-x(i,3*0+2))^2 + ...
                                                    (x(i,3*j+3)-x(i,3*0+3))^2);
                                disti = disti + 1;
                            end
                        end
                        for ja = 1:n-1
                            for jb = 1:n-1
                                if interi <= maxlength && ja ~= jb
                                    inters(interi) = sqrt((x(i,3*ja+1)-x(i,3*jb+1))^2 + ...
                                                          (x(i,3*ja+2)-x(i,3*jb+2))^2 + ...
                                                          (x(i,3*ja+3)-x(i,3*jb+3))^2);
                                    interi = interi + 1;
                                end
                            end
                        end
                    end
                end
            end
            if t <= 10 || mod(t,floor(steps/10)) == 0
                fprintf("Core=%d n=%d N0=%d dt=%d t0=%d: " + ...
                        "t=%d/%d N1=%d V1=%.20f ER=%.20f.\n", ...
                        corei,n,N0,dt,t0,t,steps,N1,V1,ER);
            end
        end
        E = mean(ERlist(floor(0.75*steps):end));
        dE = std(ERlist(floor(0.75*steps):end))/sqrt(length(ERlist(floor(0.75*steps):end)));
        runtime = floor(toc(tp));


        %print results
        fprintf("Core=%d n=%d N0=%d dt=%d t0=%d: E=%.20f±%.20e.\n",corei,n,N0,dt,t0,E,dE); %attributes
        dataPF{corei} = [dataPF{corei}; [n,N0,dt,t0,E,dE,Vsol,runtime,datai]]; %attributes
        int = fopen(main+"/intermediate.txt", 'a');
        fprintf(int, "[%d,%d,%d,%d,%.20f,%.20e,%.20f,%d];\n",n,N0,dt,t0,E,dE,Vsol,runtime); %attributes
        fseek(int,0,'eof');
        fclose(int);
        X = zeros(n,1); Y = zeros(n,1); Z = zeros(n,1);
        for j = 0:n-1
            X(j+1) = xsol(3*j+1);
            Y(j+1) = xsol(3*j+2);
            Z(j+1) = xsol(3*j+3);
        end
        F = fopen(main+"/"+sub+sprintf("/Ar%dH+.xyz",n-1), 'w');
        fprintf(F, "%d\n\nH %.20f %.20f %.20f\n", n, X(1), Y(1), Z(1));
        for j = 2:n
            fprintf(F, "Ar %.20f %.20f %.20f\n", X(j), Y(j), Z(j));
        end
        fclose(F);
        F = fopen(main+"/"+sub+"/aenergy.txt", 'w');
        fprintf(F,"%.20f\n%.20f ± %.20e\n",Vsol,E,dE);
        fclose(F);


        %plot results
        set(groot,'defaultAxesTickLabelInterpreter','latex');
        set(groot,'defaulttextinterpreter','latex');
        set(groot,'defaultLegendInterpreter','latex');
        set(groot,'defaultAxesFontSize',50); %get(groot,'factory')
        set(groot,'defaultAxesLineWidth',10);
        set(groot,'defaultLineLineWidth',10);
        set(groot,'defaultLineMarkerSize',100);
        set(groot,'defaultErrorbarLineWidth',10);
        set(groot,'defaultErrorbarMarkerSize',100);
        set(groot,'defaultAxesView',[0,90]);
        set(groot,'defaultAxesBox','on');
        set(groot,'defaultTextFontSize',50);
        set(groot,'defaultConstantlineLineWidth',10);
        set(groot,'defaultFigurePosition',[790 1 1267 1173]);
        figure; view(45,30); hold on;
        plot3(X(2:end),Y(2:end),Z(2:end),'.');
        plot3(X(1),Y(1),Z(1),'.');
        xlabel("$x$ (Bohr)"); ylabel("$y$ (Bohr)"); zlabel("$z$ (Bohr)"); 
        saveas(gcf,main+"/"+sub+"/xyz.png");
        hold off;
        figure; hold on;
        plot((1:length(ERlist))*2*t0/length(ERlist),ERlist);
        yline(E);
        xlabel("Imaginary time, $\tau$ (a.u.)");
        ylabel("Reference Energy, $E_R$ (Hartree)");
        saveas(gcf, main+"/"+sub+"/benergy.png");
        hold off;
        figure; hold on;
        dists = dists(dists~=0);
        [counts,edges] = histcounts(dists);
        counts = counts / (sum(counts)*(edges(2)-edges(1)));
        plot((edges(2:end)+edges(1:end-1))/2, counts);
        xlabel("Distance between atom and ion, $r_{ArH^+}$ (Bohr)");
        ylabel(["Probability density, $P(r_{ArH^+})$ (Bohr$^{-1}$)", ...
                "calculated over all live walkers"]);
        saveas(gcf, main+"/"+sub+"/dists.png");
        hold off;
        figure; hold on;
        inters = inters(inters~=0);
        [counts,edges] = histcounts(inters);
        counts = counts / (sum(counts)*(edges(2)-edges(1)));
        plot((edges(2:end)+edges(1:end-1))/2, counts);
        xlabel("Distance between two atoms, $r_{ArAr}$ (Bohr)");
        ylabel(["Probability density, $P(r_{ArAr})$ (Bohr$^{-1}$)", ...
                "calculated over all live walkers"]);
        saveas(gcf, main+"/"+sub+"/inters.png");
        hold off;
    end
    fprintf("Core %d finished.\n",corei);
end


%gather data
%fprintf(f,"%%n,N0,dt,t0,E,dE,Vsol,runtime\n"); %attributes
energies = cell(ln,1);
for i = 1:ln
    data(i,:) = [ns(i),N0s(i),dts(i),t0s(i),0,0,100,0]; %attributes
end
for corei = 1:cores
    for d = dataPF{corei}'
        i = d(end); %datai
        data(i,5) = data(i,5) + d(5); %attributes
        if d(7) < data(i,7)
            data(i,7) = d(7);
        end
        data(i,8) = max(data(i,8),d(8));
        energies{i} = cat(2,energies{i},d(5));
    end
end
for i = 1:ln
    data(i,5) = data(i,5) / length(energies{i});
    data(i,6) = std(energies{i}) / sqrt(length(energies{i}));
end
f = fopen(main+"/data.txt",'w');
fprintf(f,"%%n,N0,dt,t0,E,dE,Vsol,runtime\n"); %attributes
fprintf(f,"data=[...\n");
for i = 1:ln
    fprintf(f, "[%d,%d,%d,%d,%.20f,%.20e,%.20f,%d];\n", data(i,:));
end
fprintf(f,"];\n");
fclose(f);
tend = toc(tstart);
fprintf("It took %.20f seconds.\n",tend);


function f = V(x,epsilon,sigma,U,re,a1,b1,a2,b2,r1,r2,A,B,C,gamma,rs,alpha)
    n = length(x)/3;
    f = 0;
    for j = 1:n-1
        r = sqrt((x(3*j+1)-x(1))^2+(x(3*j+2)-x(2))^2+(x(3*j+3)-x(3))^2);
        if r <= r1
            if r <= re
                a = a1;
                b = b1;
            elseif r >= re
                a = a2;
                b = b2;
            end
            X = ((re/r)^a)*exp(b*(re-r));
            VArp = U * (X^2 - 2*X);
        elseif r1 <= r && r <= r2
            VArp = A / (exp(gamma*(r-rs))+1) + B*r + C;
        elseif r2 <= r
            VArp = -alpha / (2*(r^4));
        end
        f = f + VArp;
    end
    for ja = 1:n-2
        for jb = ja+1:n-1
            r = sqrt((x(3*ja+1)-x(3*jb+1))^2+(x(3*ja+2)-x(3*jb+2))^2+(x(3*ja+3)-x(3*jb+3))^2);
            f = f + 4*epsilon*((sigma/r)^12-(sigma/r)^6);
        end
    end
end