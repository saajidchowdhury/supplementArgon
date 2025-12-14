setenv('TZ', 'America/New_York');
fclose('all');
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
main = string(datetime('now','Format','M-d-y@HH.mm.ss'))+"csa2B";
mkdir(main);
tstart = tic;


%general parameters
ns = [4:13]; %[4:58] number of particles in paper, spread across different runs
cores = 10; %28 or 60 cores in paper
coresPerN = 1; %25 or more cores per cluster size in paper, spread across different runs
% parpool('Processes',cores);                            %uncomment parpool
ln = length(ns);
runs = ln*coresPerN;


%basin hopping parameters
niter = 100; %2000000 iterations in paper
R = 20;
debug = false;


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
        runtimes((i-1)*coresPerN + j) = ns(i)^2;
        jobs{(i-1)*coresPerN + j} = [ns(i),i]; %attributes
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
    fprintf(". Runtime = %d seconds.\n", totaltimes(i));
end
maxtime = max(totaltimes);
fprintf("It will take %d seconds.\n",maxtime);


%allocate and send jobs
data = zeros(ln,3); %attributes
dataPF = cell(cores,1);
parfor corei = 1:cores
    ass = assignments{corei};
    [jobn,~] = size(ass);
    for jobi = 1:jobn
        tp = tic;
        job = ass(jobi,:);
        n = job(1); datai = job(end); %attributes
        sub = sprintf("n=%dcore=%djob=%d", n, corei, jobi);
        mkdir(main+"/"+sub);

        
        %basin hopping
        fprintf("Core=%d n=%d: Starting basin.\n",corei,n);
        f = @(x)V(cat(2,xfrozen,x),epsilon,sigma,U,re,a1,b1,a2,b2,r1,r2,A,B,C,gamma,rs,alpha);
        x = zeros(1,3*(n-1)); %in Bohr
        [xsol,Vsol] = bh(f,x,R,niter,main+"/"+sub,debug);
        xsol = cat(2,xfrozen,xsol);
        fprintf("Core=%d n=%d: Finished basin. Vmin = %.20f.\n",corei,n,Vsol);
        runtime = floor(toc(tp));


        %print results
        dataPF{corei} = [dataPF{corei}; [n,Vsol,runtime,datai]]; %attributes
        int = fopen(main+"/intermediate.txt", 'a');
        fprintf(int, "[%d,%.20f,%d];\n",n,Vsol,runtime); %attributes
        fseek(int,0,'eof');
        fclose(int);
        F = fopen(main+"/"+sub+sprintf("/Ar%dH+%.6f.xyz",n-1,Vsol), 'w');
        fprintf(F, "%d\n\n", n);
        fprintf(F, "H %.20f %.20f %.20f\n", xsol(1), xsol(2), xsol(3));
        for j = 2:n
            fprintf(F, "Ar %.20f %.20f %.20f\n", xsol(3*(j-1)+1), xsol(3*(j-1)+2), xsol(3*(j-1)+3));
        end
        fclose(F);
    end
end


%gather data
%fprintf(f,"%%n,Vsol,runtime\n"); %attributes
for i = 1:ln
    data(i,:) = [ns(i),100,0]; %attributes
end
for corei = 1:cores
    for d = dataPF{corei}'
        i = d(end); %datai
        if d(2) < data(i,2)
            data(i,2) = d(2); %attributes
        end
        data(i,3) = max(data(i,3),d(3));
    end
end
f = fopen(main+"/data.txt",'w');
fprintf(f,"%%n,Vsol,runtime\n"); %attributes
fprintf(f,"data=[...\n");
for i = 1:ln
    fprintf(f, "[%d,%.20f,%d];\n", data(i,:));
end
fprintf(f,"];\n");
fclose(f);
tend = toc(tstart);
fprintf("It took %.20f seconds.\n",tend);


function [g,h] = V(x,epsilon,sigma,U,re,a1,b1,a2,b2,r1,r2,A,B,C,gamma,rs,alpha)
    g = V2B(x,epsilon,sigma,U,re,a1,b1,a2,b2,r1,r2,A,B,C,gamma,rs,alpha);
    if nargout == 2
        h = grad(x,epsilon,sigma,U,re,a1,b1,a2,b2,r1,r2,A,B,C,gamma,rs,alpha);
        h = h(4:end);
    end
end

function [xmin,Vmin] = bh(V,x,R,niter,dir,debug)
    options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','off','CheckGradients',false);
    n = length(x)/3;
    xmin = x; Vmin = 1e9;
    banksize = 50;
    bank = zeros(banksize,3*n+1);
    firstbank = zeros(banksize,3*n+1);
    totalmins = 0;
    if debug; fprintf("Starting.\n"); end
    while totalmins < niter
        for i = banksize-50+1:banksize
            still = true;
            while still
                for j = 1:n
                    x(3*(j-1)+1) = (2*rand-1) * R;
                    x(3*(j-1)+2) = (2*rand-1) * R;
                    x(3*(j-1)+3) = (2*rand-1) * R;
                    while sqrt((x(3*(j-1)+1)-x(3*(j-1)+1))^2 + ...
                            (x(3*(j-1)+2)-x(3*(j-1)+2))^2 + ...
                            (x(3*(j-1)+3)-x(3*(j-1)+3))^2) > R
                        x(3*(j-1)+1) = (2*rand-1) * R;
                        x(3*(j-1)+2) = (2*rand-1) * R;
                        x(3*(j-1)+3) = (2*rand-1) * R;
                    end
                end
                try
                    [x,Vx] = fminunc(V,x,options);
                    still = false;
                catch
                    fprintf("Linesearch error. Trying again.\n");
                end
            end
            if Vx < Vmin
                xmin = x; Vmin = Vx;
            end
            bank(i,:) = cat(2,x,Vx);
            firstbank(i,:) = cat(2,x,Vx);
        end
        Dave = 0;
        for i = 1:banksize-1
            for j = i+1:banksize
                if debug; fprintf("Bank distance is for example %d.\n",D(bank(i,1:end-1),bank(j,1:end-1))); end
                Dave = Dave + D(bank(i,1:end-1),bank(j,1:end-1));
            end
        end
        Dave = Dave / (banksize*(banksize-1)/2);
        Dcut = Dave / 2;
        if debug; fprintf("Calculated Dcut = %.2f.\n",Dcut); end
        mins = 0;
        for iter = 1:3
            flag = zeros(banksize,1);
            if debug; fprintf("Started iteration %d.\n",iter); end
            while min(flag) == 0
                counter = 0; i = 0;
                while counter < 20
                    if debug; fprintf("The counter is %d. i is %d.\n",counter,i); end
                    if min(flag) == 1
                        break;
                    end
                    i = mod(i, banksize) + 1;
                    if flag(i) == 0
                        flag(i) = 1;
                        counter = counter + 1;
                        x = bank(i,1:end-1); Vx = bank(i,end);
                        for j = 1:30
                            y = x;
                            still = true;
                            while still
                                try 
                                    if 1 <= j && j <= 20
                                        y = mate(y,firstbank,bank);
                                    elseif 21 <= j && j <= 30
                                        y = perturb(y);
                                    end
                                    [y,Vy] = fminunc(V,y,options);
                                    if Vy < Vmin
                                        xmin = y;
                                        Vmin = Vy;
                                        fprintf("Found Vmin = %.20f.\n",Vmin);
                                    end
                                    still = false;
                                catch
                                    fprintf("Linesearch error. Trying again.\n");
                                end
                            end
                            mins = mins + 1;
                            totalmins = totalmins + 1;
                            kA = 0; DalphaA = 1e9; kworst = 0; Eworst = -1e9;
                            for k = 1:banksize
                                if D(y,bank(k,1:end-1)) < DalphaA
                                    kA = k; DalphaA = D(y,bank(k,1:end-1));
                                end
                                if bank(k,end) > Eworst
                                    kworst = k; Eworst = bank(k,end);
                                end
                            end
                            if DalphaA < Dcut
                                if Vy < bank(kA,end)
                                    bank(kA,:) = cat(2,y,Vy);
                                    flag(kA) = 0;
                                end
                            elseif DalphaA >= Dcut
                                if Vy < bank(kworst,end)
                                    if debug; fprintf("Perturbed bank(%d). DalphaA = %d. Replacing bank(%d) V=%.6f with V=%.6f.\n", ...
                                    i,floor(DalphaA),kworst,bank(kworst,3*n+1),Vy); end
                                    bank(kworst,:) = cat(2,y,Vy);
                                    flag(kworst) = 0;
                                end
                            end
                        end
                    end
                end
                % if Dcut > Dave/5
                %     Dcut = Dave/2 * 0.4^(mins/10000);
                % end
            end
            fprintf("Total minimizations = %d.\n",totalmins);
        end
        banksize = banksize + 50;
    end

    sortrows(bank,3*n+1,'descend');
    G = fopen(dir+sprintf("/movie.xyz"), 'w');
    si = size(bank);
    for i = 1:min(100,si(1))
        F = fopen(dir+sprintf("/Ar%dH+bankV%.6f.xyz", n, bank(i,end)), 'w');
        fprintf(F, "%d\n\n", n+1);
        fprintf(F, "H %.20f %.20f %.20f\n", 0, 0, 0);
        for j = 1:n
            fprintf(F, "Ar %.20f %.20f %.20f\n", bank(i,3*(j-1)+1), bank(i,3*(j-1)+2), bank(i,3*(j-1)+3));
        end
        fclose(F);
        fprintf(G, "%d\n\n", n+1);
        fprintf(G, "H %.20f %.20f %.20f\n", 0, 0, 0);
        for j = 1:n
            fprintf(G, "Ar %.20f %.20f %.20f\n", bank(i,3*(j-1)+1), bank(i,3*(j-1)+2), bank(i,3*(j-1)+3));
        end
    end
    fclose(G);
end

function z = mate(y,firstbank,bank)
    n = length(y)/3;
    s = y;
    banksize = size(firstbank);
    banksize = banksize(1);
    if rand < 0.2
        k = firstbank(randi(banksize),:);
    else
        k = bank(randi(banksize),:);
    end
    k = k(1:end-1);
    nmove = floor((0.25+rand*0.25)*n);
    [qs,~] = qr(randn(3));
    [qk,~] = qr(randn(3));
    for i = 1:n
        s(3*(i-1)+1:3*(i-1)+3) = qs * s(3*(i-1)+1:3*(i-1)+3)';
        k(3*(i-1)+1:3*(i-1)+3) = qk * k(3*(i-1)+1:3*(i-1)+3)';
    end
    zs = s(3:3:end);
    zk = k(3:3:end);
    [~,inds] = sort(zs);
    [~,indk] = sort(zk);
    z = s;
    for i = 1:nmove
        z(3*(inds(i)-1)+1:3*(inds(i)-1)+3) = k(3*(indk(i)-1)+1:3*(indk(i)-1)+3);
    end
    for i = 1:n
        z(3*(inds(i)-1)+1:3*(inds(i)-1)+3) = qs' * z(3*(inds(i)-1)+1:3*(inds(i)-1)+3)';
    end
end

function z = perturb(y)
    n = length(y)/3;
    d = zeros(n,1);
    for i = 1:n
        d(i) = sqrt(y(3*(i-1)+1)^2 + y(3*(i-1)+2)^2 + y(3*(i-1)+3)^2);
    end
    [~,ind] = sort(d,'descend');
    if n >= 11
        ind = ind(1:10);
    end
    nmove = floor(rand*5) + 1;
    z = y;
    dmax = min(max(d),16);
    for i = 1:nmove
        u = rand;
        theta = acos(1-2*u);
        phi = 2*pi*rand;
        index = ind(randi(length(ind)));
        z(3*(index-1)+1) = dmax*sin(theta)*cos(phi);
        z(3*(index-1)+2) = dmax*sin(theta)*sin(phi);
        z(3*(index-1)+3) = dmax*cos(theta);
    end
end

function D = D(k,kp)
    n = length(k)/3;
    r1 = 7.55;
    r2 = 13.3;
    hk1 = zeros(n,1);
    hk2 = zeros(n,1);
    hkp1 = zeros(n,1);
    hkp2 = zeros(n,1);
    ck1 = zeros(n,1);
    ck2 = zeros(n,1);
    ckp1 = zeros(n,1);
    ckp2 = zeros(n,1);
    for i = 1:n-1
        for j = i+1:n
            r = sqrt((k(3*(i-1)+1)-k(3*(j-1)+1))^2 + ...
                (k(3*(i-1)+2)-k(3*(j-1)+2))^2 + ...
                (k(3*(i-1)+3)-k(3*(j-1)+3))^2);
            if 0 < r && r < r1
                ck1(i) = ck1(i) + 1;
                ck1(j) = ck1(j) + 1;
            elseif r1 < r && r < r2
                ck2(i) = ck2(i) + 1;
                ck2(j) = ck2(j) + 1;
            end
            r = sqrt((kp(3*(i-1)+1)-kp(3*(j-1)+1))^2 + ...
                (kp(3*(i-1)+2)-kp(3*(j-1)+2))^2 + ...
                (kp(3*(i-1)+3)-kp(3*(j-1)+3))^2);
            if 0 < r && r < r1
                ckp1(i) = ckp1(i) + 1;
                ckp1(j) = ckp1(j) + 1;
            elseif r1 < r && r < r2
                ckp2(i) = ckp2(i) + 1;
                ckp2(j) = ckp2(j) + 1;
            end
        end
    end
    for i = 1:n
        hk1(ck1(i)+1) = hk1(ck1(i)+1) + 1;
        hk2(ck2(i)+1) = hk2(ck2(i)+1) + 1;
        hkp1(ckp1(i)+1) = hkp1(ckp1(i)+1) + 1;
        hkp2(ckp2(i)+1) = hkp2(ckp2(i)+1) + 1;
    end
    D = 0;
    for i = 0:n-1
        D = D + i * (2*abs(hk1(i+1)-hkp1(i+1)) + abs(hk2(i+1)-hkp2(i+1)));
    end
end

function f = V2B(x,epsilon,sigma,U,re,a1,b1,a2,b2,r1,r2,A,B,C,gamma,rs,alpha)
    threshold = 20;
    n = length(x)/3;
    f = 0;
    for j = 2:n
        for k = 1:3
            if threshold < x(3*(j-1)+k)
                f = f + 100*(x(3*(j-1)+k)-threshold);
            elseif x(3*(j-1)+k) < -threshold
                f = f + 100*(-threshold-x(3*(j-1)+k));
            end
        end
    end
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

function g = grad(x,epsilon,sigma,U,re,a1,b1,a2,b2,r1,r2,A,B,C,gamma,rs,alpha)
    threshold = 20;
    n = length(x)/3;
    g = zeros(3*n,1);
    for j = 2:n
        for k = 1:3
            if threshold < x(3*(j-1)+k)
                g(3*(j-1)+k) = g(3*(j-1)+k) + 100;
            elseif x(3*(j-1)+k) < -threshold
                g(3*(j-1)+k) = g(3*(j-1)+k) - 100;
            end
        end
    end
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
            dXdr = a*((re/r)^(a-1))*(-re/(r^2))*exp(b*(re-r)) + ((re/r)^a)*(-b)*exp(b*(re-r));
            dVdr = U * (2*X - 2) * dXdr;
        elseif r1 <= r && r <= r2
            dVdr = - A / (exp(gamma*(r-rs)) + 1)^2 * gamma*exp(gamma*(r-rs)) + B;
        elseif r2 <= r
            dVdr = 4 * alpha / (2*(r^5));
        end
        for k = 0:2
            g(3*j+k+1) = g(3*j+k+1) + dVdr * (x(3*j+k+1)-x(k+1)) / r;
            g(k+1) = g(k+1) + dVdr * (x(k+1)-x(3*j+k+1)) / r;
        end
    end
    for ja = 1:n-2
        for jb = ja+1:n-1
            r = sqrt((x(3*ja+1)-x(3*jb+1))^2+(x(3*ja+2)-x(3*jb+2))^2+(x(3*ja+3)-x(3*jb+3))^2);
            dVdr = 4*epsilon*(12*((sigma/r)^11)*(-sigma/(r^2))-6*((sigma/r)^5)*(-sigma/(r^2)));
            for k = 0:2
                g(3*ja+k+1) = g(3*ja+k+1) + dVdr * (x(3*ja+k+1)-x(3*jb+k+1)) / r;
                g(3*jb+k+1) = g(3*jb+k+1) + dVdr * (x(3*jb+k+1)-x(3*ja+k+1)) / r;
            end
        end
    end
end
