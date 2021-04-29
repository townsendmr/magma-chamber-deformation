function [Q] = heat_conduction_chamber(maxn,a,c,dr,kappa,rho,cp,Tb)


global storeTime storeTemp storeSumk storeSumk_2 storeSumk_old storeSumk_2_old lengthTime maxTime
global storeSumk_oldold storeSumk_2_oldold switch_Tprofile

% geometry
r     = a+dr:dr:a+dr; % radial coordinate (m)


% time grid
time = storeTime;
time_index= length(time);

if maxTime < time(end)
    maxTime = time(end);
    storeSumk_oldold=storeSumk_old;
    storeSumk_2_oldold=storeSumk_2_old;
    storeSumk_old=storeSumk;
    storeSumk_2_old=storeSumk_2;
end
% initial and boundary conditions
Tr0   = storeTemp;
% if maxTime > time(end)
%     storeSumk_old=storeSumk;
%     storeSumk_2_old=storeSumk_2;
% end
%if time_index>2
%     display([num2str(time_index)  ' ; ' num2str(time(end-2)) ' ; ' num2str(time(end-1)) ' ; ' num2str(time(end))]);
%     display(['index: ' num2str(time_index)]);%  ' ; ' num2str(time(end-2)) ' ; ' num2str(time(end-1)) ' ; ' num2str(time(end))]);

%end
% pick a time
current_time    = time(end);

% WHAT I SHOULD DO HERE IS DO ALL THE STEP EXCEPT THE LAST 2 SO THAT IF
% ADAPTIVE TIME STEP I AM COVERED
if time_index == 2
    % first sum over n
    
    sumn = 0;
    for n=1:maxn
        % sum over k within first sum over n
        sumk=0;
        for k=1:time_index-1
            past_time       = time(k);
            past_time_delta = time(k+1);
            
            mk       = (Tr0(k+1) - Tr0(k))/(time(k+1)-time(k));
            pk       = Tr0(k) - mk*time(k);
            
            termk = mk*c^4/(kappa^2*n^4*pi^4) ...
                * (exp(-kappa*(n*pi/c)^2*(current_time-past_time))*(1-kappa*(n*pi/c)^2*past_time) ...
                + exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta))*(kappa*(n*pi/c)^2*past_time_delta-1)) ...
                + pk*c^2/(kappa*n^2*pi^2) ...
                * (exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) ...
                - exp(-kappa*(n*pi/c)^2*(current_time-past_time)));
            sumk = sumk + termk;
        end
        lengthTime=2;
        storeSumk(n)=sumk;
        termn = n.*sin(n.*pi.*(r-a)./c).*sumk;
        sumn  = sumn + termn;
    end
    
    
    term1 = -4.*pi.*kappa.*a./(2.*r.*c^2).*(-1).*sumn;
    
    
    % second sum over n
    sumn = 0;
    for n=1:maxn
        % sum over k within second sum over n
        sumk = 0;
        for k=1:time_index-1
            past_time       = time(k);
            past_time_delta = time(k+1);
            termk = exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) ...
                - exp(-kappa*(n*pi/c)^2*(current_time-past_time));
            sumk = sumk + termk;
        end
        storeSumk_2(n)=sumk;
        termn = 1./n.*sin(n.*pi.*(r-a)./c).*cos(n*pi).*sumk;
        sumn  = sumn + termn;
    end
    
    
    term2 = -4.*pi.*kappa.*(a+c)./(2.*kappa.*pi^2.*r).*Tb.*(1).*sumn;
elseif time_index > 2 && time_index>lengthTime
    % first sum over n
  %  display(['index: ' num2str(time_index) ' ; lengthTime : ' num2str(lengthTime) ' ; time-dt : ' num2str(time(time_index-1)) ' ; time : ' num2str(time(time_index)) ]);%  ' ; ' num2str(time(end-2)) ' ; ' num2str(time(end-1)) ' ; ' num2str(time(end))]);
    
    sumn = 0;
    %display([num2str(time_index)  ' ; ' num2str(time(end-2)) ' ; ' num2str(time(end-1)) ' ; ' num2str(time(end))]);
    for n=1:maxn
        % sum over k within first sum over n
        
        temp=storeSumk(n)*exp(-kappa*(n*pi/c)^2*(current_time-time(time_index-1)));
        
        past_time       = time(time_index-1);
        past_time_delta = time(time_index);
        
        mk       = (Tr0(time_index) - Tr0(time_index-1))/(time(time_index)-time(time_index-1));
        pk       = Tr0(time_index-1) - mk*time(time_index-1);
        
        termk = mk*c^4/(kappa^2*n^4*pi^4) ...
            * (exp(-kappa*(n*pi/c)^2*(current_time-past_time))*(1-kappa*(n*pi/c)^2*past_time) ...
            + exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta))*(kappa*(n*pi/c)^2*past_time_delta-1)) ...
            + pk*c^2/(kappa*n^2*pi^2) ...
            * (exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) ...
            - exp(-kappa*(n*pi/c)^2*(current_time-past_time)));
        sumk = temp + termk;
        %storeSumk_old(n)=storeSumk(n);
        storeSumk(n)=sumk;
        termn = n.*sin(n.*pi.*(r-a)./c).*sumk;
        sumn  = sumn + termn;
    end
    
    
    term1 = -4.*pi.*kappa.*a./(2.*r.*c^2).*(-1).*sumn;
    
   % display(['CH Term 1 = ' num2str(term1)]);

    % second sum over n
    sumn = 0;
    for n=1:maxn
        temp=storeSumk_2(n)*exp(-kappa*(n*pi/c)^2*(current_time-time(time_index-1)));
        past_time       = time(time_index-1);
        past_time_delta = time(time_index);
        termk = exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) ...
            - exp(-kappa*(n*pi/c)^2*(current_time-past_time));
        sumk = temp + termk;
        %storeSumk_2_old(n)=storeSumk_2(n);
        storeSumk_2(n)=sumk;
        termn = 1./n.*sin(n.*pi.*(r-a)./c).*cos(n*pi).*sumk;
        sumn  = sumn + termn;
    end
    term2 = -4.*pi.*kappa.*(a+c)./(2.*kappa.*pi^2.*r).*Tb.*(1).*sumn;
     %   display(['CH Term 2 = ' num2str(term2)]);
        switch_Tprofile=0;
elseif time_index > 2 && time_index==lengthTime
    sumn = 0;
 %   display([num2str(time_index)  ' ; ' num2str(time(end-2)) ' ; ' num2str(time(end-1)) ' ; ' num2str(time(end))]);
 %   display('TROUBLE');
 %   display(['index: ' num2str(time_index) ' ; lengthTime : ' num2str(lengthTime) ' ; time-dt : ' num2str(time(end-1)) ' ; time : ' num2str(time(end)) ]);%  ' ; ' num2str(time(end-2)) ' ; ' num2str(time(end-1)) ' ; ' num2str(time(end))]);
 %   display(['SwithTProfile = ' num2str(switch_Tprofile)])
    sumn = 0;
    %display([num2str(time_index)  ' ; ' num2str(time(end-2)) ' ; ' num2str(time(end-1)) ' ; ' num2str(time(end))]);
    for n=1:maxn
        % sum over k within first sum over n
        if switch_Tprofile==0
            temp=storeSumk_old(n)*exp(-kappa*(n*pi/c)^2*(current_time-time(end-1)));
            past_time       = time(time_index-1);
            past_time_delta = time(time_index);
            
            mk       = (Tr0(time_index) - Tr0(time_index-1))/(time(time_index)-time(time_index-1));
            pk       = Tr0(time_index-1) - mk*time(time_index-1);
            
            termk = mk*c^4/(kappa^2*n^4*pi^4) ...
                * (exp(-kappa*(n*pi/c)^2*(current_time-past_time))*(1-kappa*(n*pi/c)^2*past_time) ...
                + exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta))*(kappa*(n*pi/c)^2*past_time_delta-1)) ...
                + pk*c^2/(kappa*n^2*pi^2) ...
                * (exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) ...
                - exp(-kappa*(n*pi/c)^2*(current_time-past_time)));
            sumk = temp + termk;
            %storeSumk_old(n)=storeSumk(n);
            storeSumk(n)=sumk;
        else
            sumk=storeSumk_old(n);
            %display('I am here.........')
            %switch_Tprofile=0;
        end
        
        termn = n.*sin(n.*pi.*(r-a)./c).*sumk;
        sumn  = sumn + termn;
    end
    
    
    term1 = -4.*pi.*kappa.*a./(2.*r.*c^2).*(-1).*sumn;
    %display(['CH Term 1 = ' num2str(term1)]);
    
    % second sum over n
    sumn = 0;
    for n=1:maxn
        if switch_Tprofile==0
            
            temp=storeSumk_2_old(n)*exp(-kappa*(n*pi/c)^2*(current_time-time(end-1)));
            past_time       = time(end-1);
            past_time_delta = time(end);
            termk = exp(-kappa*(n*pi/c)^2*(current_time-past_time_delta)) ...
                - exp(-kappa*(n*pi/c)^2*(current_time-past_time));
            sumk = temp + termk;
            %storeSumk_2_old(n)=storeSumk_2(n);
            storeSumk_2(n)=sumk;
            
        else
            sumk=storeSumk_2_old(n);
            %switch_Tprofile=0;
        end
        termn = 1./n.*sin(n.*pi.*(r-a)./c).*cos(n*pi).*sumk;
        sumn  = sumn + termn;
    end
    term2 = -4.*pi.*kappa.*(a+c)./(2.*kappa.*pi^2.*r).*Tb.*(1).*sumn;
   % display(['CH Term 2 = ' num2str(term2)]);

elseif time_index > 2 && time_index<lengthTime
    sumn = 0;
    switch_Tprofile=1;
   % display([num2str(time_index)  ' ; ' num2str(time(end-2)) ' ; ' num2str(time(end-1)) ' ; ' num2str(time(end))]);
   % display(['TROUBLE' num2str(switch_Tprofile)]);
   % display(['index: ' num2str(time_index) ' ; lengthTime : ' num2str(lengthTime) ' ; time-dt : ' num2str(time(end-1)) ' ; time : ' num2str(time(end)) ]);%  ' ; ' num2str(time(end-2)) ' ; ' num2str(time(end-1)) ' ; ' num2str(time(end))]);
    
    sumn = 0;
    %display([num2str(time_index)  ' ; ' num2str(time(end-2)) ' ; ' num2str(time(end-1)) ' ; ' num2str(time(end))]);
    for n=1:maxn
        
        termn = n.*sin(n.*pi.*(r-a)./c).*storeSumk_old(n);
        %storeSumk(n)=storeSumk_old(n);
        sumn  = sumn + termn;
    end
    
    
    term1 = -4.*pi.*kappa.*a./(2.*r.*c^2).*(-1).*sumn;
    
   % display(['CH Term 1 = ' num2str(term1)]);

    % second sum over n
    sumn = 0;
    for n=1:maxn
        
        termn = 1./n.*sin(n.*pi.*(r-a)./c).*cos(n*pi).*storeSumk_2_old(n);
        %storeSumk_2(n)=storeSumk_2_old(n);
        sumn  = sumn + termn;
    end
    term2 = -4.*pi.*kappa.*(a+c)./(2.*kappa.*pi^2.*r).*Tb.*(1).*sumn;
    %     storeSumk_old=storeSumk;
    %     storeSumk_2_old=storeSumk_2;
    %    display(['CH Term 2 = ' num2str(term2)]);

elseif time_index <2
    term1 =zeros(size(r));
    term2 =zeros(size(r));
end

% third sum over n
Ta = Tr0(time_index);
sumn = 0;
for n=1:maxn
%     termn =  sin(n.*pi.*(r-a)./c).*exp(-kappa.*(n.*pi./c)^2.*current_time) ...
%         * ((a*(a+c)*Ta - a*(a+c)*Tb)/(n*pi)*(1-cos(n*pi)) ... %SIGN!!
%         + ((a+c)*Tb -a*Ta)/(n*pi)*(-(a+c)*cos(n*pi) + a));
    termn =  sin(n.*pi.*(r-a)./c).*exp(-kappa.*(n.*pi./c)^2.*current_time) ...
        * ((a*(a+c)*Ta - a*(a+c)*Tb)/(n*pi)*(1-cos(n*pi)) ... %SIGN!!
        + ((a+c)*Tb -a*Ta)/(n*pi)*(-(a+c)*cos(n*pi) + a));
    
    sumn = sumn + termn;
end
term3 = 2./(r.*c).*sumn;
%display(['CH Term 3 = ' num2str(term3)]);

Trt = term1 + term2 + term3;
%display(['Trt CH = ' num2str(Trt)]);

lengthTime=time_index;
small_q = -kappa*rho*cp*(Trt-Tr0(time_index))/dr;
surface_area_chamber=4*pi*a^2;
Q = small_q*surface_area_chamber;
