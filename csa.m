load("interpolation.mat");
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
main = string(datetime('now','Format','M-d-y@HH.mm.ss'))+"csa";
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
R = 10;
debug = false;


%potential parameters
a0 = 0.529177210544; %Angstroms
xfrozen = [1.5056,0,0,0,0,0,-1.5056,0,0] / a0; %Bohr
EhK = 3.1577502480398e5; %Kelvin
Ktocm = 0.695034800486127; %cm^{-1}


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
        f = @(x)V(cat(2,xfrozen*a0,x),rlist,thetalist,V_2dlist,potararlist);
        x = zeros(1,3*(n-3)); %in Angstrom
        [xsol,Vsol] = bh(f,x,R,niter,main+"/"+sub,debug);
        xsol = xsol / a0; %now, in Bohr
        xsol = cat(2,xfrozen,xsol);
        Vsol = Vsol / Ktocm / EhK;
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
        fprintf(F, "Ar %.20f %.20f %.20f\n", xsol(1), xsol(2), xsol(3));
        fprintf(F, "H %.20f %.20f %.20f\n", xsol(4), xsol(5), xsol(6));
        fprintf(F, "Ar %.20f %.20f %.20f\n", xsol(7), xsol(8), xsol(9));
        for j = 4:n
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


function [g,h] = V(x,rlist,thetalist,V_2dlist,potararlist)
    g = VJuditlinear(x,rlist,thetalist,V_2dlist,potararlist);
    if nargout == 2
        h = gradVJudit(x,rlist,thetalist,V_2dlist,potararlist);
        h = h(10:end);
    end
end

function [xmin,Vmin] = bh(V,x,R,niter,dir,debug)
    a0 = 0.529177210544; %Angstroms
    EhK = 3.1577502480398e5; %Kelvin
    Ktocm = 0.695034800486127; %cm^{-1}
    options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','off');
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
                    if debug; fprintf("Linesearch error. Trying again.\n"); end
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
                                        if debug; fprintf("Found Vmin = %.20f.\n",Vmin/EhK/Ktocm); end
                                    end
                                    still = false;
                                catch
                                    if debug; fprintf("Linesearch error. Trying again.\n"); end
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
                                    i,floor(DalphaA),kworst,bank(kworst,3*n+1)/EhK/Ktocm,Vy/EhK/Ktocm); end
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
                % if debug; fprintf("Dcut = %.2f.\n",Dcut); end
            end
            fprintf("Total minimizations = %d.\n",totalmins);
        end
        banksize = banksize + 50;
        fprintf("Total minimizations = %d.\n",totalmins);
    end

    % if debug
        sortrows(bank,3*n+1);
        G = fopen(dir+sprintf("/movie.xyz"), 'w');
        for i = 1:size(bank,1)
            F = fopen(dir+sprintf("/Ar%dH+bankV%.6f.xyz", n+2, bank(i,end)/EhK/Ktocm), 'w');
            fprintf(F, "%d\n\n", n+3);
            fprintf(F, "Ar %.20f %.20f %.20f\n", 1.5056/a0, 0, 0);
            fprintf(F, "H %.20f %.20f %.20f\n", 0, 0, 0);
            fprintf(F, "Ar %.20f %.20f %.20f\n", -1.5056/a0, 0, 0);
            for j = 1:n
                fprintf(F, "Ar %.20f %.20f %.20f\n", bank(i,3*(j-1)+1)/a0, bank(i,3*(j-1)+2)/a0, bank(i,3*(j-1)+3)/a0);
            end
            fclose(F);
            fprintf(G, "%d\n\n", n+3);
            fprintf(G, "Ar %.20f %.20f %.20f\n", 1.5056/a0, 0, 0);
            fprintf(G, "H %.20f %.20f %.20f\n", 0, 0, 0);
            fprintf(G, "Ar %.20f %.20f %.20f\n", -1.5056/a0, 0, 0);
            for j = 1:n
                fprintf(G, "Ar %.20f %.20f %.20f\n", bank(i,3*(j-1)+1)/a0, bank(i,3*(j-1)+2)/a0, bank(i,3*(j-1)+3)/a0);
            end
        end
        fclose(G);
    % end
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

function z = perturb(y) %change this to 1-5 atoms at a time 
    % n = length(y)/3;
    % r1 = 4;
    % c = zeros(n,1);
    % for i = 1:n-1
    %     for j = i+1:n
    %         r = sqrt((y(3*(i-1)+1)-y(3*(j-1)+1))^2 + ...
    %             (y(3*(i-1)+2)-y(3*(j-1)+2))^2 + ...
    %             (y(3*(i-1)+3)-y(3*(j-1)+3))^2);
    %         if r1 <= r
    %             c(i) = c(i) + 1;
    %             c(j) = c(j) + 1;
    %         end
    %     end
    % end
    % [csorted,ind] = sort(c);
    % indlowest = ind(csorted==min(csorted));
    % if min(csorted) ~= max(csorted)
    %     indsecondlowest = ind(csorted<=min(csorted(csorted>min(csorted))));
    % else
    %     indsecondlowest = indlowest;
    % end
    % ilowest = indlowest(randi(length(indlowest)));
    % isecondlowest = indsecondlowest(randi(length(indsecondlowest)));
    % z = y;
    % z(3*(ilowest-1)+1:3*(ilowest-1)+3) = z(3*(isecondlowest-1)+1:3*(isecondlowest-1)+3) + (2*rand-1) * r1;
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
    dmax = min(max(d),8.3);
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
    r1 = 4;
    r2 = 7;
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
    rmsxk = 0;
    rmsyzk = 0;
    rmsxkp = 0;
    rmsyzkp = 0;
    for i = 1:n
        rmsxk = rmsxk + min(k(3*(i-1)+1)^2, 8.3^2);
        rmsyzk = rmsyzk + min(k(3*(i-1)+2)^2 + k(3*(i-1)+3)^2, 6.4^2);
        rmsxkp = rmsxkp + min(kp(3*(i-1)+1)^2, 8.3^2);
        rmsyzkp = rmsyzkp + min(kp(3*(i-1)+2)^2 + kp(3*(i-1)+3)^2, 6.4^2);
    end
    rmsxk = sqrt(rmsxk);
    rmsyzk = sqrt(rmsyzk);
    rmsxkp = sqrt(rmsxkp);
    rmsyzkp = sqrt(rmsyzkp);
    % D = D + 50*sqrt((rmsxk-rmsxkp)^2+(rmsyzk-rmsyzkp)^2);
    D = D + 1000*abs(rmsyzk-rmsyzkp);
end

%potential (in cm^-1, given coordinates in Angstroms)
function totenergy = VJuditlinear(coordinate,rlist,thetalist,V_2dlist,potararlist)
    n = length(coordinate) / 3;
    totenergy = 0;
    vector1 = zeros(3,1);
    vector1(1) = coordinate(1) - coordinate(7);
    vector1(2) = coordinate(2) - coordinate(8);
    vector1(3) = coordinate(3) - coordinate(9);
    a = sqrt(vector1(1)^2 + vector1(2)^2 + vector1(3)^2);
    for i = 4:n
        vector2(1) = coordinate(4) - coordinate(3*(i-1)+1);
        vector2(2) = coordinate(5) - coordinate(3*(i-1)+2);
        vector2(3) = coordinate(6) - coordinate(3*(i-1)+3);
        r = sqrt(vector2(1)^2 + vector2(2)^2 + vector2(3)^2);
        cosang = (vector1(1)*vector2(1) + vector1(2)*vector2(2) + vector1(3)*vector2(3)) / (a*r);
        v = V_2dlinear(r,cosang, rlist,thetalist,V_2dlist);
        totenergy = totenergy + v;
    end
    for i = 4:n-1
        for j = i+1:n
            r = sqrt((coordinate(3*(i-1)+1)-coordinate(3*(j-1)+1))^2 + ...
                (coordinate(3*(i-1)+2)-coordinate(3*(j-1)+2))^2 + ...
                (coordinate(3*(i-1)+3)-coordinate(3*(j-1)+3))^2);
            v = potararlinear(r, rlist,potararlist);
            totenergy = totenergy + v;
        end
    end
    totenergy = totenergy + 38155.94046159532119*(n-4);
end

function v = V_2dlinear(r,cost,rlist,thetalist,V_2dlist)
    nr = length(rlist);
    nt = length(thetalist);
    theta = acos(cost);
    j = floor((theta-thetalist(1))/(thetalist(end)-thetalist(1))*(nt-1)) + 1;
    assert(1 <= j && j <= nt-1, "j out of range");
    mt = (theta-thetalist(j))/(thetalist(j+1)-thetalist(j));
    if rlist(1) < r && r < rlist(end)
        i = floor((r-rlist(1))/(rlist(end)-rlist(1))*(nr-1)) + 1;
        assert(1 <= i && i <= nr-1, "i out of range");
        mr = (r-rlist(i))/(rlist(i+1)-rlist(i));
        v = V_2dlist(i,j) + mr*(V_2dlist(i+1,j)-V_2dlist(i,j)) + ...
                            mt*(V_2dlist(i,j+1)-V_2dlist(i,j));
    else
        assert(rlist(end) <= r, "r is tiny");
        v = V_2dlist(nr,j) + (r-rlist(nr)) * 10.0;
    end
end

function v = potararlinear(r,rlist,potararlist)
    nr = length(rlist);
    if rlist(1) < r && r < rlist(end)
        i = floor((r-rlist(1))/(rlist(end)-rlist(1))*(nr-1)) + 1;
        assert(1 <= i && i <= nr-1, "i out of range");
        mr = (r-rlist(i))/(rlist(i+1)-rlist(i));
        v = potararlist(i) + mr*(potararlist(i+1)-potararlist(i));
    else
        assert(rlist(end) <= r, "r is tiny");
        v = potararlist(nr) + (r-rlist(nr)) * 1.0;
    end
end

function grad = gradVJudit(coordinate,rlist,thetalist,V_2dlist,potararlist)
    n = length(coordinate) / 3;
    grad = zeros(3*n,1);
    vector1 = zeros(3,1);
    vector1(1) = coordinate(1) - coordinate(7);
    vector1(2) = coordinate(2) - coordinate(8);
    vector1(3) = coordinate(3) - coordinate(9);
    a = sqrt(vector1(1)^2 + vector1(2)^2 + vector1(3)^2);
    for i = 4:n
        vector2(1) = coordinate(4) - coordinate(3*(i-1)+1);
        vector2(2) = coordinate(5) - coordinate(3*(i-1)+2);
        vector2(3) = coordinate(6) - coordinate(3*(i-1)+3);
        r = sqrt(vector2(1)^2+vector2(2)^2+vector2(3)^2);
        cost = (vector1(1)*vector2(1) + ...
               vector1(2)*vector2(2) + ...
               vector1(3)*vector2(3)) / (a*r);
        it = floor((acos(cost)-thetalist(1))/(thetalist(end)-thetalist(1))*(length(thetalist)-1)) + 1;
        if rlist(1) < r && r < rlist(end)
            ir = floor((r-rlist(1))/(rlist(end)-rlist(1))*(length(rlist)-1)) + 1;
            dV_2ddr = (V_2dlist(ir+1,it) - V_2dlist(ir,it)) / (rlist(ir+1)-rlist(ir));
            dV_2ddc = (V_2dlist(ir,it+1)-V_2dlist(ir,it)) / (cos(thetalist(it+1))-cos(thetalist(it)));
        else
            dV_2ddr = 10.0;
            dV_2ddc = 0.0;
        end
        grad(3*(i-1)+1) = grad(3*(i-1)+1) + dV_2ddr*coordinate(3*(i-1)+1)/r + dV_2ddc*(-1/r+coordinate(3*(i-1)+1)^2/r^3);
        grad(3*(i-1)+2) = grad(3*(i-1)+2) + dV_2ddr*coordinate(3*(i-1)+2)/r + dV_2ddc*(coordinate(3*(i-1)+1)*coordinate(3*(i-1)+2)/r^3);
        grad(3*(i-1)+3) = grad(3*(i-1)+3) + dV_2ddr*coordinate(3*(i-1)+3)/r + dV_2ddc*(coordinate(3*(i-1)+1)*coordinate(3*(i-1)+3)/r^3);
    end
    for i = 4:n
        for j = 4:n
            if i ~= j
                r = sqrt((coordinate(3*(i-1)+1)-coordinate(3*(j-1)+1))^2 + ...
                    (coordinate(3*(i-1)+2)-coordinate(3*(j-1)+2))^2 + ...
                    (coordinate(3*(i-1)+3)-coordinate(3*(j-1)+3))^2);
                if rlist(1) < r && r < rlist(end)
                    ir = floor((r-rlist(1))/(rlist(end)-rlist(1))*(length(rlist)-1)) + 1;
                    dpotarardr = (potararlist(ir+1)-potararlist(ir)) / (rlist(ir+1)-rlist(ir));
                else
                    dpotarardr = 1.0;
                end
                grad(3*(i-1)+1) = grad(3*(i-1)+1) + dpotarardr*(coordinate(3*(i-1)+1)-coordinate(3*(j-1)+1))/r;
                grad(3*(i-1)+2) = grad(3*(i-1)+2) + dpotarardr*(coordinate(3*(i-1)+2)-coordinate(3*(j-1)+2))/r;
                grad(3*(i-1)+3) = grad(3*(i-1)+3) + dpotarardr*(coordinate(3*(i-1)+3)-coordinate(3*(j-1)+3))/r;
            end
        end
    end
end
