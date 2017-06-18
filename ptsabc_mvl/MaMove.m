function [msp] = MaMove(fn,mn,fsp,msp,femass,mamass,d,lb,ub,pm)
%MAMOVE Summary of this function goes here
%   Detailed explanation goes here
    %Preliminaries
    dt=zeros(1,mn);
    % Scale for distance
    scale=(-lb(1)+ub(1));%/2;
    [Indb, pvkd] = find(mamass>=median(mamass));  % Male spiders above median
    for i=1:mn
        if ismember(i,Indb)     % Spider above the median
            % Start looking for a female with stronger vibration
            for q1=1:fn
                if femass(q1)>mamass(i)
                    % Calculate the distance
                    dt(q1)=norm(msp(i,:)-fsp(q1,:));
                else
                    dt(q1)=0;
                end
            end
            % Choose the shortest distance
            [pvkd1,Ind,val] = find(dt);   % Choose where the distance in non zero
            [pvkd2,Imin] = min(val);      % Get the shortest distance
            Ish = Ind(Imin);
            % Update moves
            if isempty(val)
                Vib=0;
                spaux=zeros(1,d);
            else
                dt=dt./scale;
                Vib = 2*femass(Ish)*exp(-(rand*dt(Ish).^2));
                spaux=fsp(Ish,:);
            end
            delta = 2*rand(1,d)-.5;
            tmpf = 2*pm.*(rand(1,d)-0.5);
            msp(i,:) = msp(i,:)+Vib*(spaux-msp(i,:)).*delta+tmpf;
        else % de aqui para abajo falta
            %% Spider below median, go to weigthed mean
            % Generate the weighted mean
            spdpos = [fsp' msp']';
            spdwei = [femass' mamass']';
            weigth = repmat(spdwei,1,d);
            dim = find(size(spdpos)~=1,1);
            wmean = sum(weigth.*spdpos,dim)./sum(weigth,dim);
            %% Move
            delta = 2*rand(1,d)-.5;
            tmpf = 2*pm.*(rand(1,d)-0.5);
            msp(i,:) = msp(i,:)+(wmean-msp(i,:)).*delta+tmpf;
        end
    end
    % Check limits
%         for i=1:d
%             for q1=1:mn
%                 if msp(q1,i)<lb(i), msp(q1,i)=lb(i)+(ub(i)-lb(i)).*rand(1,1); end
%                 if msp(q1,i)==lb(i), msp(q1,i)=lb(i); end
% 
%                 if msp(q1,i)>ub(i), msp(q1,i)=lb(i)+(ub(i)-lb(i)).*rand(1,1); end
%                 if msp(q1,i)==ub(i), msp(q1,i)=ub(i); end
%             end
%         end
for i= 1:d
          for q1=1:mn
     if  (exp(j*pi/4)<=msp(q1,i)<exp(3*j*pi/4))                                           %(s(m)>=exp(j*pi/4) )&& (s(m)<exp(3*j*pi/4))
        msp(q1,i)=j;
      else if msp(q1,i)>=exp(3*j*pi/4)&& msp(q1,i)<exp(5*j*pi/4)
         msp(q1,i)=-1;
          else if msp(q1,i)>=exp(5*j*pi/4)&& msp(q1,i)<exp(7*j*pi/4)
                msp(q1,i)=-j;
              else
                  msp(q1,i)=1;
             end
          end
     end  
          end
      end
      







end