%2D Brownian motion via particle collisions. Collisions handled as
%described in Wang (1992) or Foerster 1994. Derek W. Harrison August 25, 2015.

clear 
clc
tic

%Parameters
ma=2.0;                               %Mass of particle 1 (the large particle)
mpi=0.9;                            %Mass of particles other than particle 1
Ra=0.7;                             %Radius of particle 1 (m)
Rpi=0.2;                            %Radii of particle i (m)
Dwally= 3;                          %Size of system in positive and negative y direction (m)
Dwallx = 3;                         %Size of system in positive and negative x direction (m)
Nboxya = 2*Ra;                      %Range in the y direction (positive and negative) (m) relative to particle 1 in which possible particle collisions are checked. If particle j, which may collide with 1, is not within range Nboxya of 1 collision time is not checked.
Nboxxa = 2*Ra;                      %Range in the x direction (positive and negative) (m) relative to particle 1 in which possible particle collisions are checked. If particle j, which may collide with 1, is not within range Nboxxa of 1 collision time is not checked.
Nboxyi = 0.5;                       %Range in the y direction (positive and negative) (m) relative to particle i in which possible particle collisions are checked. If particle j, which may collide with i, is not within range Nboxyi of i collision time is not checked.
Nboxxi = 0.5;                       %Range in the x direction (positive and negative) (m) relative to particle i in which possible particle collisions are checked. If particle j, which may collide with i, is not within range Nboxxi of i collision time is not checked.
Vixmax = 2;                       %Maximum initial x velocity component for particles (m) other than particle 1
Viymax = 2;                       %Maximum initial y velocity component for particles (m) other than particle 1

Nts = 20000;                        %Number of timesteps. Note: this value is not equal to the actual number of timesteps taken in simulation due to the fact that collision times are smaller than dt. 
Nfs = 5000;                         %Number of frames for video
Ns = 15;                             %Number of simulations. Note: Setting value above to one leads to large simulations times. In order to save computation time (and also prevent videos from exceeding size limits) the parameters writevideo and Plot are set to 0.
to = 0;                             %Initial time (s)
tf = 10;                            %Final time (s)

dt = (tf-to)/Nts;                   %Timestepsize (s)
dtfs = (tf-to)/Nfs;                 %Plots are stored at intervals dtfs

e  = 1.0;                           %normal restitution coefficient dogzs
mu = 0.00;                          %Friction coefficient
Bo = 1.0;                           %Tangential restitution

error=1e-25;                        %Small negligible number used in the denominator of various equations to prevent equations from blowing up
smallerr = 1e-5;                    %Number used for preventing particles from sticking with walls and from particles going through eachother (in case of identical collision times for various particles)
maxfractionp=0.75;                  %Maximum area available for particles other than particle 1 (the large particle)
fracfactor = 0.7;                   %Fraction factor. It is multiplied with maxfractionp to vary amount of particles placed in system. Note: High fractions lead to large Monte Carlo simulation times for generating initial particle positions
writevideo = 0;                     %Button for video. Set to 1 to generate video of simulation. Set to any other arbitrary number to prevent generation of video.
Plot = 0;                           %Generate plot of simulation at everytimestep
frameps = 120;                       %Set framerate (fps) for video
%End adjustable parameters

A = 4*Dwally*Dwallx;                %Area of system box
frac = maxfractionp*fracfactor;     %Fraction of available area for particles
A1 = A - 4*pi*Ra^2;                 %Area of particles available for all particles but particle 1 (the large particle)

Ia=2/5*ma*Ra^2;                     %Moment of inertia of particle 1
Ib = 2/5*mpi*Rpi^2;                 %Moment of inertia of all particles but particle 1

Npi = floor(A1*frac/(pi*Rpi^2));    %Number of particles placed in system other than particle 1 (large particle)
Np = Npi+1;                         %Total number of particles in system

Xo = zeros(Np,1);
Yo = zeros(Np,1);

Vox = zeros(Np,1);
Voy = zeros(Np,1);
Wox = zeros(Np,1);
Woy = zeros(Np,1);

Ri = zeros(Np,1);
mi = zeros(Np,1);
Ii = zeros(Np,1);

%Generating Initial conditions for large particle
Xo(1,1) = 0;
Yo(1,1) = 0;
Vox(1,1) = 0;
Voy(1,1) = 0;
Wox(1,1) = 0;
Woy(1,1) = 0;

%Initializing radius, mass and moment of inertia vectors all particles but
%particle 1
Ri(1) = Ra;mi(1) = ma;Ii(1) = Ia;
for i=2:Np
    Ri(i) = Rpi;
    mi(i) = mpi;
    Ii(i) = Ib;    
end

%Multisimulations start here
dRaveragesum = zeros(Nts+5,1);
timestepmin = zeros(Ns,1);

if writevideo==1
    writerObj = VideoWriter('C:\Users\d-w-h\Desktop\Home\Sim videos.avi','Motion JPEG AVI');
    writerObj.FrameRate = frameps;
    open(writerObj);
end

%Multisimulation loop starts here
for simc = 1:Ns
    %Generating initial velocities of all particles but particle 1
    for i=2:Np
        Vox(i)=(2*rand-1)*Vixmax;
        Voy(i)=(2*rand-1)*Viymax;
    end

    %Generating Initial positions via Monte Carlo simulation for all particles
    %but particles 1
    for i=2:Np 
        Xo(i) = (rand*2-1)*(Dwallx-Rpi*1.1);
        Yo(i) = (rand*2-1)*(Dwally-Rpi*1.1);
        there_is_overlap = true;
        while there_is_overlap == true
            there_is_overlap = false;
            for j=1:i-1            
                dXvec = sqrt((Xo(i)-Xo(j))^2);
                dYvec = sqrt((Yo(i)-Yo(j))^2);
                if dXvec^2+dYvec^2 < 1.01*(Ri(i)+Ri(j))^2
                    there_is_overlap = true;
                end
            end
            if there_is_overlap == true
                Xo(i) = (rand*2-1)*(Dwallx-Rpi*1.1);
                Yo(i) = (rand*2-1)*(Dwally-Rpi*1.1);
            end
        end
    end

    %Main simulation code starts here
    if Plot == 1
        plot(Xo,Yo,'b.','MarkerSize',60)
        hold on
        plot(Xo(1),Yo(1),'r.','MarkerSize',60)
        hold off
        title('Brownian motion');
        minY = -Dwally;
        maxY = Dwally;
        minX = -Dwallx;
        maxX = Dwallx;
        xlabel('x-coordinate')
        ylabel('y-coordinate')
        set(gca,'Ylim',[minY maxY])
        set(gca,'Xlim',[minX maxX])
    end

    %Some parameters
    X = Xo;
    Y = Yo;
    Vx = Vox;
    Vy = Voy;
    Wx = Wox;
    Wy = Woy;
    rwallx = 2*Np;              %Rightwall lable
    lwallx = 2*Np+1;            %Leftwall lable
    uwally = 2*Np+2;            %Upperwall lable
    lwally = 2*Np+3;            %lowerwall lable

    %The simulation iterations begin here
    time = 0;
    timestep = 1;
    CheckingColltime = zeros(Nts*10,1);
    Colltimesarray = zeros(Nts*10,1);
    doubletimesarray = zeros(Nts*10,1);
    timearray = zeros(10*Nts,1);
    dtarray = zeros(10*Nts,1);
    collpart1counter = 0;
    framecounter = 1;
    Overlaparray = zeros(Nts*10,1);
    X1array = zeros(Nts*10,1);
    Y1array = zeros(Nts*10,1);
    X1array(1) = X(1);
    Y1array(1) = Y(1);

    while time < tf
        colltime = dt;
        collisionwithwall = 0;
        collisionwithparticle = 0;
        colltimecounter = 0;
        %Checking collision times and finding collision partner     
        for i=1:Np
            for j = i+1:Np
                tab = dt;
                rab =[X(i)-X(j);Y(i)-Y(j)];vab=[Vx(i)-Vx(j);Vy(i)-Vy(j)];vab2=vab'*vab;
                if (i > 1) && (abs(X(i)-X(j)) < Nboxxi) && (abs(Y(i)-Y(j)) < Nboxyi)
                    Disc = (rab'*vab)^2-vab2*(rab'*rab-(Ri(i)+Ri(j))^2);
                    %Check overlap
                    if sqrt(rab'*rab) < (Ri(i)+Ri(j))
                        Overlaparray(timestep) = Overlaparray(timestep)+1;
                    end
                    if (Disc > 0) 
                        tab = (-rab'*vab-sqrt(Disc))/vab2;
                    end
                    if (tab <= colltime) && (tab >= 0)
                        colltime = tab;
                        PartnerA = i;PartnerB=j;
                        collisionwithparticle=1;
                        collisionwithwall=0;
                        colltimecounter = colltimecounter+1;
                        CheckingColltime(colltimecounter)=colltime;
                    end
                elseif (i == 1) && (abs(X(i)-X(j)) < Nboxxa) && (abs(Y(i)-Y(j)) < Nboxya)
                    Disc = (rab'*vab)^2-vab2*(rab'*rab-(Ri(i)+Ri(j))^2);
                    if sqrt(rab'*rab) < (Ri(i)+Ri(j))
                        Overlaparray(timestep) = Overlaparray(timestep)+1;
                    end
                    if (Disc > 0) 
                        tab = (-rab'*vab-sqrt(Disc))/vab2;
                    end
                    if (tab <= colltime) && (tab >= 0)
                        colltime = tab;
                        PartnerA = i;PartnerB=j;
                        collisionwithparticle=1;
                        collisionwithwall=0;
                        colltimecounter = colltimecounter+1;
                        CheckingColltime(colltimecounter)=colltime;
                    end
                end
            end

            colltimerwallx = (Dwallx-Ri(i)-X(i))/(Vx(i)+error);
            if (colltimerwallx >= 0) && (colltimerwallx <= colltime)
                colltime = colltimerwallx;
                PartnerA = i;
                PartnerB = rwallx;
                collisionwithwall = 1;
                collisionwithparticle = 0;
                CheckingColltime(i)=colltime;
                colltimecounter = colltimecounter+1;
                CheckingColltime(colltimecounter)=colltime;        
            end

            colltimelwallx = -(X(i)-Ri(i)+Dwallx)/(Vx(i)+error);
            if (colltimelwallx >= 0) && (colltimelwallx <= colltime)
                colltime = colltimelwallx;
                PartnerA = i;
                PartnerB = lwallx;
                collisionwithwall = 1;
                collisionwithparticle = 0;
                CheckingColltime(i)=colltime;
                colltimecounter = colltimecounter+1;
                CheckingColltime(colltimecounter)=colltime;
            end

            colltimeuwally = (Dwally-Ri(i)-Y(i))/(Vy(i)+error);
            if (colltimeuwally >= 0) && (colltimeuwally <= colltime)
                colltime = colltimeuwally;
                PartnerA = i;
                PartnerB = uwally;
                collisionwithwall = 1;
                collisionwithparticle = 0;
                CheckingColltime(i)=colltime;
                colltimecounter = colltimecounter+1;
                CheckingColltime(colltimecounter)=colltime;
            end

            colltimelwally = -(Y(i)-Ri(i)+Dwally)/(Vy(i)+error);
            if (colltimelwally >= 0) && (colltimelwally <= colltime)
                colltime = colltimelwally;
                PartnerA = i;
                PartnerB = lwally;
                collisionwithwall = 1;
                collisionwithparticle = 0;
                CheckingColltime(i)=colltime;
                colltimecounter = colltimecounter+1;
                CheckingColltime(colltimecounter)=colltime;
            end   
        end

        if (PartnerA == 1) && (PartnerB <= Np)
            collpart1counter = collpart1counter+1;
        end

        doubletimes = 0;
        CheckingColltime(CheckingColltime == 0) = [];
        sizeCheck = size(CheckingColltime);
        for i = 1:sizeCheck(1)
            for j = i+1:sizeCheck(1)
                if (CheckingColltime(i) == CheckingColltime(j)) && (CheckingColltime(i) < dt)
                    doubletimes=doubletimes+1;
                end
            end    
        end

        doubletimesarray(timestep) = doubletimes;
        Colltimesarray(timestep) = colltime;

        %Update positions and velocities 
        if colltime < dt 
            %Update positions to the point of collision minus some small, negligible,  numerical value to
            %prevent particles from getting stuck in a wall
            X = X + Vx*colltime*(1-smallerr);
            Y = Y + Vy*colltime*(1-smallerr);

            if (collisionwithparticle == 1)        
                %Update velocities of colliding particles
                ra = [X(PartnerA);Y(PartnerA)];va = [Vx(PartnerA);Vy(PartnerA)];
                rb = [X(PartnerB);Y(PartnerB)];vb = [Vx(PartnerB);Vy(PartnerB)];              
                n = (ra-rb)/sqrt((ra-rb)'*(ra-rb));
                vab=[Vx(PartnerA)-Vx(PartnerB);Vy(PartnerA)-Vy(PartnerB)];
                wa = [Wx(PartnerA);Wy(PartnerA)];
                wb = [Wx(PartnerB);Wy(PartnerB)];           
                RiWi = Ri(PartnerA)*wa + Ri(PartnerB)*wb;
                crossRiWin = cross([RiWi;0], [n;0]);
                vab = vab - crossRiWin(1:2);

                B1 = 7/2*(1/mi(PartnerA)+1/mi(PartnerB));
                B2 = 1/mi(PartnerA)+1/mi(PartnerB);    

                t = (vab - n*(vab'*n))/sqrt((vab - n*(vab'*n))'*(vab - n*(vab'*n))+error);        
                Jn = -(1+e)*(vab'*n)/B2;

                stickyslide = (1+Bo)*(vab'*t)/Jn/B1;

                if mu < stickyslide %Sliding 
                    Jt = -mu*Jn;                
                elseif mu >= stickyslide %Sticking
                    Jt = -(1+Bo)*(vab'*t)/B1;      
                end

                J = Jn*n+Jt*t;
                Vx(PartnerA) = J(1)/mi(PartnerA)+Vx(PartnerA);Vy(PartnerA) = J(2)/mi(PartnerA)+Vy(PartnerA);
                Vx(PartnerB) = -J(1)/mi(PartnerB)+Vx(PartnerB);Vy(PartnerB) = -J(2)/mi(PartnerB)+Vy(PartnerB);

                crossnJ = cross([n;0],[-J;0]);

                Wx(PartnerA) = -crossnJ(1)*Ri(PartnerA)/Ii(PartnerA)+Wx(PartnerA);
                Wy(PartnerA) = -crossnJ(2)*Ri(PartnerA)/Ii(PartnerA)+Wy(PartnerA);
                Wx(PartnerB) = -crossnJ(1)*Ri(PartnerB)/Ii(PartnerB)+Wx(PartnerB);
                Wy(PartnerB) = -crossnJ(2)*Ri(PartnerB)/Ii(PartnerB)+Wy(PartnerB);      

            elseif (collisionwithwall == 1)
                if (PartnerB == rwallx) || (PartnerB == lwallx)
                    %Update velocity of particle A
                    Vx(PartnerA) = -Vx(PartnerA);
                end
                if (PartnerB == uwally) || (PartnerB == lwally)
                    %Update velocity of particle A
                    Vy(PartnerA) = -Vy(PartnerA);
                end     
            end      
            time = time + colltime;
            dtarray(timestep) = colltime;
            timestep=timestep+1;
            timearray(timestep) = time;
            X1array(timestep) = X(1);
            Y1array(timestep) = Y(1);
        elseif colltime >= dt %Update positions and velocities using dt
            X = X + Vx*dt;
            Y = Y + Vy*dt;
            time = time + dt
            dtarray(timestep) = dt;
            timestep=timestep+1;
            timearray(timestep) = time;
            X1array(timestep) = X(1);
            Y1array(timestep) = Y(1); 
            simc
        end

        dummy = floor(time/dtfs)+1; %For storing plots at intervals dtfs
        if (dummy == framecounter) && (Plot == 1) && (simc == 1)  
            plot(X,Y,'b.','MarkerSize',60)
            hold on
            plot(X1array(1:timestep),Y1array(1:timestep),X(1),Y(1),'r.','MarkerSize',175)
            hold off
            title('Brownian motion');
            minY = -Dwally;
            maxY = Dwally;
            minX = -Dwallx;
            maxX = Dwallx;
            xlabel('x-coordinate')
            ylabel('y-coordinate')
            set(gca,'Ylim',[minY maxY])
            set(gca,'Xlim',[minX maxX])
            set(gca,'nextplot','replacechildren');
            set(gcf,'Renderer','zbuffer')
            %set(gcf,'units','normalized','outerposition',[0.2 0.10 0.56 1.6*0.56])  
            frame = getframe(gcf); 
            framecounter=framecounter+1;
            if writevideo == 1      
                writeVideo(writerObj,frame);
            end
        end
        
    end%End of single simulation loop

    X1array = X1array(1:timestep);
    Y1array = Y1array(1:timestep);
    timearray = timearray(1:timestep);
    Darray = zeros(size(timearray));
    dRsquare = zeros(size(timearray));
    DarrayL = zeros(size(timearray));
    collpart1counter;
    timestepmin(simc) = timestep;

    if writevideo == 1
        close(writerObj);
    end

    %Calculating diffusion coefficient
    for i = 2:timestep
        dRsquare(i,1) = (X1array(i)-X1array(1))^2+(Y1array(i)-Y1array(1))^2;
        Darray(i,1) = ((X1array(i)-X1array(1))^2+(Y1array(i)-Y1array(1))^2)/4/timearray(i);
    end

    %Interpolation
    timearrayint = to:(tf-to)/Nts:tf;
    timearrayint = timearrayint';
    dRsquareint = interp1(timearray,dRsquare,timearrayint);

    dRaveragesum(1:Nts+1) = dRaveragesum(1:Nts+1) + dRsquareint;
end%End of Main multisimulation loop

dRaverage = dRaveragesum/simc;
Darrayavg = zeros(size(timearrayint));
dRaverage = dRaverage(1:Nts+1);

%Calculating average diffusion coefficient
for i = 2:Nts+1
    Darrayavg(i,1) = dRaverage(i)/3/timearrayint(i);
end

%Plotting some data
Dmean = mean(Darray);
Dmeanavg = mean(Darrayavg);
figure
plot(timearray,dRsquare,timearray,4*Dmean*timearray,timearrayint,dRaverage, timearrayint,4*Dmeanavg*timearrayint)

figure
plot(timearrayint,dRaverage, timearrayint,4*Dmeanavg*timearrayint)

%Export Data to excel
A = [X1array,Y1array,timearray,Darray];
B = [timearrayint,dRaverage,Darrayavg];
filename = 'C:\Users\d-w-h\Desktop\Home\Brownian Motion Simulation results\BrownianMotionSim36Data.xlsx';
xlswrite(filename,A,1,'A1:D25000')
xlswrite(filename,B,2,'A1:C25000')

toc  