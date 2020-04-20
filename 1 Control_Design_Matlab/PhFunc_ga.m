function y = PhFunc_ga(PPP)
global PP BB Kh A1

% 2 inputs
P1=PPP(1,1);
P2=PPP(1,2);
P3=PPP(1,3);
P4=PPP(1,4);
Phat = [P1 P2;P3 P4];

n = sqrt(numel(Phat));

% eigs = eig(Phat);
% penalti1 = zeros(1,n);
% penalti2 = zeros(1,n);
% penalti3 = zeros(1,n);
% for i = 1:n
%     if  isreal(eigs(i)) == 0
%         penalti1(i) = 1;
%     elseif real(eigs(i)) < 0
%         penalti2(i) = 1;
%     end
% end

K = Phat\Kh;
eigs1 = (eig(A1+BB*K));
for i = 1:4
    if  real(eigs1(i)) >= 0
        penalti3(i) = 100;
%     elseif real(eigs1(i)) < -100
%         penalti3(i) = 100;
    else
        penalti3(i) = 0;
    end
end
% penalti = sum(penalti1)+sum(penalti2)+sum(penalti3);
penalti = sum(penalti3);
y = max(max(abs((PP*BB-BB*Phat))))+penalti;