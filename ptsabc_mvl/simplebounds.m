function s=simplebounds(s,lb,ub)
ns_tmp=s;
        I=ns_tmp<lb;
        ns_tmp(I)=lb(I);
        I1=ns_tmp>ub;
        ns_tmp(I1)=ub(I1);
        s=ns_tmp;
        for m= 1:4
     if  (exp(j*pi/4)<=s(m)<exp(3*j*pi/4))                                           %(s(m)>=exp(j*pi/4) )&& (s(m)<exp(3*j*pi/4))
         s(m)=j;
     else if s(m)>=exp(3*j*pi/4)&& s(m)<exp(5*j*pi/4)
         s(m)=-1;
         else if s(m)>=exp(5*j*pi/4)&& s(m)<exp(7*j*pi/4)
               s(m)=-j;
             else
                 s(m)=1;
             end
         end
     end            
        end