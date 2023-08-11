%Rate independent crystal plasticity Taylor code (Full Constraints)
%Darshan Chalapathi, IIT Madras
%email-id: subbi1darshan@gmail.com
%%
clear
clc
% fp1 = fopen('VP.txt','w');
fp2 = fopen('final_tex.txt','w');

%Read data from the Taylor.in file
%read_data() is the function file and must be present in the same folder
[eulerangles,iopt,nslip,ns1,ns2,bs,CRSS,nsteps,tdot,L,ihard,lin_hard_const,tau0_1,tau1_1,theta0_1,theta1_1,tau0_2,tau1_2,theta0_2,theta1_2] = read_data();
if iopt == 1
    sname = 'fcc-';
elseif iopt == 2
    sname = 'bcc-2-';
elseif iopt == 3
    sname = 'bcc-1-';
end
ngrain = size(eulerangles,1);

%Normalising the planes and directions
temp1 = sqrt(bs(1,1)^2 + bs(1,2)^2 + bs(1,3)^2);
bs = bs./temp1;
if iopt == 1 || iopt == 3
    temp2 = sqrt(ns1(1,1)^2 + ns1(1,2)^2 + ns1(1,3)^2);
    ns1 = ns1./temp2;
    CRSS(1:ngrain,1:nslip) = CRSS;
else
    temp2 = sqrt(ns1(1,1)^2 + ns1(1,2)^2 + ns1(1,3)^2);
    ns1 = ns1./temp2;
    temp3 = sqrt(ns2(1,1)^2 + ns2(1,2)^2 + ns2(1,3)^2);
    ns2 = ns2./temp3;
    CRSS_1 = CRSS(1);
    CRSS_2 = CRSS(2);
    CRSS(1:ngrain,1:nslip/2) = CRSS_1;
    CRSS(1:ngrain,(nslip/2)+1:nslip) = CRSS_2;
end

%Determining the Schmid tensor
schmidtensor(1:3,1:3,1:nslip) = 0;
if iopt == 1 || iopt == 3
    for s = 1:nslip/2
        for i = 1:3
            for j = 1:3
                schmidtensor(i,j,s) = (bs(s,i)*ns1(s,j));
                schmidtensor(i,j,(s+12)) = -schmidtensor(i,j,s);
            end
        end
    end
else
    for s = 1:nslip/4
        for i = 1:3
            for j = 1:3
                schmidtensor(i,j,s) = (bs(s,i)*ns1(s,j));
                schmidtensor(i,j,(s+12)) = -(bs(s,i)*ns1(s,j));
                schmidtensor(i,j,(s+24)) = (bs(s,i)*ns2(s,j));
                schmidtensor(i,j,(s+36)) = -(bs(s,i)*ns2(s,j));
            end
        end
    end
end

%Defining the basis tensors for vectorisation process
basis(:,:,1) = [-1 0 0;0 -1 0;0 0 2]./sqrt(6);
basis(:,:,2) = [-1 0 0;0 1 0;0 0 0]./sqrt(2);
basis(:,:,3) = [0 0 0;0 0 1;0 1 0]./sqrt(2);
basis(:,:,4) = [0 0 1;0 0 0;1 0 0]./sqrt(2);
basis(:,:,5) = [0 1 0;1 0 0;0 0 0]./sqrt(2);

%Determining the symmetric and anti-symmetric part of velocity gradient tensor
for i = 1:3
    for j = 1:3
        D(i,j) = 0.5*(L(i,j)+L(j,i));
        W(i,j) = 0.5*(L(i,j)-L(j,i));
    end
end

%Vectorizing the strain rate
VD(1:5) = 0;
VD(1) = ((D(3,3) - D(1,1))+(D(3,3) - D(2,2)))/sqrt(6);
VD(2) = (D(2,2)-D(1,1))/sqrt(2);
VD(3) = sqrt(2)*D(2,3);
VD(4) = sqrt(2)*D(1,3);
VD(5) = sqrt(2)*D(1,2);

acc_shear(ngrain,1) = 0;
CRSS(1:ngrain,1:nslip) = CRSS;

for istep = 1:nsteps   %Loop over all increments
    for igrain = 1:ngrain   %Loop over all orientations
        phi1 = eulerangles(igrain,1)*pi/180;
        phi = eulerangles(igrain,2)*pi/180;
        phi2 = eulerangles(igrain,3)*pi/180;
        
        %Calculation of orientation matrix
        g_phi1 = [cos(phi1) sin(phi1) 0;-sin(phi1) cos(phi1) 0;0 0 1];
        g_phi = [1 0 0;0 cos(phi) sin(phi);0 -sin(phi) cos(phi)];
        g_phi2 = [cos(phi2) sin(phi2) 0;-sin(phi2) cos(phi2) 0;0 0 1];
        rot_mat = g_phi2*g_phi*g_phi1;

        %Rotating the schmid tensors onto the sample frame by considering the
        %transpose of rotation matrix determined above.
        rot_mat = rot_mat.';
        rot_schmidtensor(1:3,1:3,1:nslip) = 0;
        
        for islip=1:nslip   %Loop over sll slip systems
            for i = 1:3
                for j = 1:3
                    for k = 1:3
                        for l = 1:3
                            rot_schmidtensor(i,j,islip) = rot_schmidtensor(i,j,islip) +(rot_mat(i,k)*rot_mat(j,l)*schmidtensor(k,l,islip));
                        end
                    end
                end
            end
            
            %Converting to symmetric and anti-symmetric Schmid tensors
            for i = 1:3
                for j = 1:3
                    P(i,j,islip) = 0.5*(rot_schmidtensor(i,j,islip)+rot_schmidtensor(j,i,islip));
                    Q(i,j,islip) = 0.5*(rot_schmidtensor(i,j,islip)-rot_schmidtensor(j,i,islip));
                end
            end
        end     %End over all slip systems   

        %Vectorizing the symmetric schmid tensor
        for i = 1:nslip
            VP(1,i) = ((P(3,3,i) - P(1,1,i))+(P(3,3,i) - P(2,2,i)))/sqrt(6);
            VP(2,i) = (P(2,2,i)-P(1,1,i))/sqrt(2);
            VP(3,i) = sqrt(2)*P(2,3,i);
            VP(4,i) = sqrt(2)*P(1,3,i);
            VP(5,i) = sqrt(2)*P(1,2,i);
            %fprintf(fp1,'%f %f %f %f %f\n',VP(1,i),VP(2,i),VP(3,i),VP(4,i),VP(5,i));
        end

        %Using simplex method
        Aeq = VP;         %Equality constraint 1
        beq = VD;         %Equality constraint 2
        lb = zeros(nslip,1);
        A = [];                     %For Ax <= B
        b = [];                     %NO inequality constraints
        for i = 1:nslip
            f(i) = CRSS(igrain,i); %Objective function - Minimisation of sum of shear rates
        end
        options = optimoptions('linprog','Algorithm','dual-simplex');
        [x,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb,[],[],'options');
        stress = lambda.eqlin;      %Stress obtained directly from lambda multipliers
        stress.' * VP;              %Projecting stress on Schmid tensor to check the validity
        stress1(igrain,:) = stress.';
        
        %Obtaining the stress tensor from vectorised tensor
        stress_tensor(1:3,1:3,igrain) = 0;
        for i = 1:5
            stress_tensor(1:3,1:3,igrain) = stress_tensor(1:3,1:3,igrain) + stress1(igrain,i)*basis(1:3,1:3,i);
        end
        
        shearrates(igrain,:) = x;
        TF(igrain)=sum(shearrates(igrain,:));

        %Neglecting very small shear rates and assigning slip system to be active or not
        n_size=size(shearrates);
        n_rows=n_size(1,1);
        n_cols=n_size(1,2);
        for p=1:n_rows
            for q=1:n_cols
                if shearrates(p,q) <= 1e-7
                    shearrates(p,q) = 0;
                    tau(p,q) = 0;       %Not active
                else
                    tau(p,q) = 1;       %Active
                end
            end
        end
        
        %Reorient the grain after deformation
        aux(1:3,1:3) = 0;
        for islip = 1:nslip
            for i = 1:3
                for j = 1:3
                    aux(i,j) = aux(i,j) + (Q(i,j,islip) * shearrates(igrain,islip));
                end
            end
        end
	
        for i = 1:3
            for j = 1:3
                omega(i,j) = (W(i,j) - aux(i,j))*tdot;
            end
        end
        
        instant_shear(igrain,:) = abs(shearrates(igrain,:).*tdot);
        delta_shear(igrain,1) = sum(abs(shearrates(igrain,:).*tdot));
        acc_shear(igrain,1) = acc_shear(igrain) + sum(abs(shearrates(igrain,:).*tdot));
        
        v(1) = omega(3,2);
        v(2) = omega(1,3);
        v(3) = omega(2,1);
        snorm = sqrt((v(1)*v(1)) + (v(2)*v(2)) + (v(3)*v(3)));
        snorm1 = tan(snorm/2);
        if snorm <= 1e-6
            snorm = 1;
        end
        for i = 1:3
            vbar(i) = snorm1 * v(i)/snorm;
        end
        snorm = (vbar(1)*vbar(1)) + (vbar(2)*vbar(2)) + (vbar(3)*vbar(3));
        th(3,2) = vbar(1);
        th(1,3) = vbar(2);
        th(2,1) = vbar(3);
        th(2,3) = -vbar(1);
        th(3,1) = -vbar(2);
        th(1,2) = -vbar(3);
        for i = 1:3
            th(i,i) = 0;
        end

        for i = 1:3
            for j = 1:3
                th2(i,j) = 0;
                for k = 1:3
                    th2(i,j) = th2(i,j) + th(i,k)*th(k,j);
                end
            end
        end
        
        iden33 = eye(3);
        for i = 1:3
            for j = 1:3
                rot(i,j) = iden33(i,j) + 2.*(th(i,j)+th2(i,j))/(1 + snorm);
            end
        end
        
        for i = 1:3
            for j = 1:3
                anew(i,j) = 0;
                for k = 1:3
                    anew(i,j) = anew(i,j) + rot(i,k)*rot_mat(k,j);
                end
            end
        end
        anew = anew.';
        
        th1 = acos(anew(3,3));
        if abs(anew(3,3) >= 0.9999)
            tm = 0;
            ph = atan2(anew(1,2),anew(1,1));
        else
            sth = sin(th1);
            tm = atan2(anew(1,3)/sth,anew(2,3)/sth);
            ph = atan2(anew(3,1)/sth,-anew(3,2)/sth);
        end
        phi1_1 = ph*180/pi;
        phi_1 = th1*180/pi;
        phi2_1 = tm*180/pi;
        eulerangles(igrain,1) = phi1_1;
        eulerangles(igrain,2) = phi_1;
        eulerangles(igrain,3) = phi2_1;
        if istep == nsteps       %Writing the reoriented Euler angles at final step
            fprintf(fp2,'%f\t%f\t%f\n',eulerangles(igrain,1),eulerangles(igrain,2),eulerangles(igrain,3));
        end
        
        %Linear hardening
        if ihard == 1
            for i = 1:nslip
                CRSS(igrain,i) = CRSS(igrain,i) + lin_hard_const * instant_shear(igrain,i);
            end
        end
        
        %Voce hardening
        hlatex = 1.0;
        if ihard == 2
           if iopt ~= 2
               for p = 1:nslip
                   for q = 1:nslip
                       if (p == q)
                          hard(p,q) = 1;
                       else
                          hard(p,q) = hlatex;
                       end
                   end
               end
               for p = 1:nslip
                   dtau = 0;
                   for q = 1:nslip
                       dtau = dtau + hard(p,q) * instant_shear(igrain,q);
                   end
                   tiny = 1e-4 * tau0_1;
                   voce = 0;
                   if (abs(theta0_1) > tiny)
                       voce = theta1_1 * delta_shear(igrain);
                       if (abs(tau1_1) > tiny)
                           fact = abs(theta0_1/tau1_1);
                           exp_ini = exp(-acc_shear(igrain)*fact);
                           exp_del = exp(-delta_shear(igrain)*fact);
                           voce = voce - (((fact*tau1_1 - theta1_1)/fact)*exp_ini*(exp_del - 1)) - ...
                                  ((theta1_1/fact)*exp_ini*(exp_del*(fact*(acc_shear(igrain)+delta_shear(igrain))+1)...
                                  -((acc_shear(igrain)*fact)+1)));
                       end
                   end
                   CRSS(igrain,p) = CRSS(igrain,p) + dtau*voce/delta_shear(igrain);
               end
            else
               for p = 1:nslip
                   for q = 1:nslip
                       if (p == q)
                          hard(p,q) = 1;
                       else
                          hard(p,q) = hlatex;
                       end
                   end
               end
               for p = 1:nslip/2
                   dtau = 0;
                   for q = 1:nslip/2
                       dtau = dtau + hard(p,q) * instant_shear(igrain,q);
                   end
                   tiny = 1e-4 * tau0_1;
                   voce = 0;
                   if (abs(theta0_1) > tiny)
                       voce = theta1_1 * delta_shear(igrain);
                       if (abs(tau1_1) > tiny)
                           fact = abs(theta0_1/tau1_1);
                           exp_ini = exp(-acc_shear(igrain)*fact);
                           exp_del = exp(-delta_shear(igrain)*fact);
                           voce = voce - (((fact*tau1_1 - theta1_1)/fact)*exp_ini*(exp_del - 1)) - ...
                                  ((theta1_1/fact)*exp_ini*(exp_del*(fact*(acc_shear(igrain)+delta_shear(igrain))+1)...
                                  -((acc_shear(igrain)*fact)+1)));
                       end
                   end
                   CRSS(igrain,p) = CRSS(igrain,p) + dtau*voce/delta_shear(igrain);
               end
               for p = (nslip/2)+1:nslip
                   dtau = 0;
                   for q = (nslip/2)+1:nslip
                       dtau = dtau + hard(p,q) * instant_shear(igrain,q);
                   end
                   tiny = 1e-4 * tau0_2;
                   voce = 0;
                   if (abs(theta0_2) > tiny)
                       voce = theta1_2 * delta_shear(igrain);
                       if (abs(tau1_2) > tiny)
                           fact = abs(theta0_2/tau1_2);
                           exp_ini = exp(-acc_shear(igrain)*fact);
                           exp_del = exp(-delta_shear(igrain)*fact);
                           voce = voce - (((fact*tau1_2 - theta1_2)/fact)*exp_ini*(exp_del - 1)) - ...
                                  ((theta1_2/fact)*exp_ini*(exp_del*(fact*(acc_shear(igrain)+delta_shear(igrain))+1)...
                                  -((acc_shear(igrain)*fact)+1)));
                       end
                   end
                   CRSS(igrain,p) = CRSS(igrain,p) + dtau*voce/delta_shear(igrain);
               end
            end   
        end
        
        %Calculating von Mises equivalent stress    
        stress_vm(igrain,istep) = sqrt(0.5*((stress_tensor(1,1,igrain)-stress_tensor(2,2,igrain))^2+...
                                  (stress_tensor(2,2,igrain)-stress_tensor(3,3,igrain))^2+(stress_tensor(3,3,igrain)-stress_tensor(1,1,igrain))^2+...
                                  6*(stress_tensor(1,2,igrain)^2+stress_tensor(2,3,igrain)^2+stress_tensor(3,1,igrain)^2)));
    end
%     if istep == 13
%         sname1 = '20';
%         sname2 = strcat(sname,sname1,'.out');
%         save(sname2,'eulerangles','-ascii')
%         sname3 = strcat(sname,sname1,'.mat');
%         save(sname3)
%     elseif istep == 30
%         sname1 = '40';
%         sname2 = strcat(sname,sname1,'.out');
%         save(sname2,'eulerangles','-ascii')
%         sname3 = strcat(sname,sname1,'.mat');
%         save(sname3)
%     elseif istep == 53
%         sname1 = '60';
%         sname2 = strcat(sname,sname1,'.out');
%         save(sname2,'eulerangles','-ascii')
%         sname3 = strcat(sname,sname1,'.mat');
%         save(sname3)
%     elseif istep == 93
%         sname1 = '80';
%         sname2 = strcat(sname,sname1,'.out');
%         save(sname2,'eulerangles','-ascii')
%         sname3 = strcat(sname,sname1,'.mat');
%         save(sname3)        
%     end
end

% fclose(fp1);
% fclose(fp2);

% Calculating Taylor factor
L_norm = 0;
for i = 1:3
    for j = 1:3
        L_sym(i,j) = 0.5*(L(i,j) + L(j,i));
        L_norm = L_norm + (L_sym(i,j)*L_sym(i,j));
    end
end
L_norm = sqrt((2/3) * L_norm);
taylor_factor = TF./L_norm;
taylor_factor = taylor_factor.';
save('taylor_factor.txt','taylor_factor','-ascii')
save('stress_values.txt','stress1','-ascii')