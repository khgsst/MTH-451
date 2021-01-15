clear, close all
dims=zeros((9+2*15)*12,1);
rhos = dims;

ii = 0;
d = 0;
for f=1:9
    m=f*10^d;
    for j = 1:12
        ii = ii+1;
        A = rand(m)*2-1;
        [~,U,P]=lu(A);
        rho = max(max(abs(U)))/max(max(abs(A)));
        rhos(ii) = rho;
        dims(ii) = m;
    end
end

for d=1:2
    for f=[1 1.1 1.2 1.5 1.7 2 2.5 3 3.5 4 5 6 7 8 9]
        m=fix(f*10^d);
        for j = 1:12
            ii = ii+1;
            A = rand(m)*2-1;
            [~,U,P]=lu(A);
            rho = max(max(abs(U)))/max(max(A));
            rhos(ii)=rho;
            dims(ii) = m;
        end
    end
end

loglog (dims,rhos,'.',1:1000, sqrt(1:1000), '--', ...
    1:1000,(1:1000).^(2/3),'-.','linewidth',2);
legend('measured','m^{1/2}','m^{2/3}');
xlabel('m')
ylabel('growth factor \rho')
title('Growth factor version dimension')
disp('hit any key to continue')
pause
close
disp('GENERATING PDF vs. \rho:')
%%
clear
figure(5)
hold on
rhos=zeros(10^6,1);
for m=[8,16,32]
    for ii = 1:10^6
        A = rand(m)*2-1;
        [~,U,P]=lu(A);
        rho = max(max(abs(U)))/max(max(abs(A)));
        rhos(ii) = rho;
    end
    x=0:0.01:10;
    N=histogram(rhos,x);
    pdf = N/(10^6*0.01);
    semilogy(x,pdf) 
    %hold on
end
grid on
xlabel('\rho')
ylabel('pdf')
legend('m=8','m=16','m=32')
disp('hit any key to continue')
pause
hold off
close
disp('Now generating 4 plots of entries of inverse L')
%%
A = rand(128)*2-1;
[L,u,p]=lu(A);
[q,~]=qr(L);
figure(1)
spy(max(abs(inv(L))));
xlabel(['max abs entry of L^{-1} = ',num2str(max(max(abs(inv(L)))))])
figure(2)
spy(abs(q)>1/sqrt(128))
xlabel('random A')
Lp = L.*sign(randn(128));
[qp,r]=qr(Lp);
figure(3)
spy(abs(inv(Lp))>1)
xlabel(['max abs entry of L^{-1} = ',num2str(max(max(abs(inv(Lp)))))])
figure(4)
spy(abs(qp)>1/sqrt(128));
xlabel('L with random sign')


