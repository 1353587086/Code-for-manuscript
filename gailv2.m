A=zeros(721,2001);
for i=1:321
    for j=1:2001
        xxx=0.25*(j-1);yyy=0.25*(i-1);
        xigema11=Kfind(diyi,xxx,yyy);xigema33=Kfind(disan3,xxx,yyy);%第一层
       if xigema33>0
        miuc=c(i,j)*1e6;miufai=fai(i,j)/180;
        Z=(miuc*cos(miufai)+sin(miufai)*(xigema11+xigema33)/2)/((xigema11-xigema33)/2)-1;
        ZZ=(((2*cos(miufai)/(xigema11-xigema33))^2)*(xigemac(i,j)^2)+((xigema11+xigema33)/(((xigema11-xigema33))*cos(miufai)-2*miuc*sin(miufai)/(xigema11-xigema33))^2)*(xigemafai(i,j)^2))^0.5;
        baita=Z/ZZ;
        Pf=1-normcdf(baita,0,1);
        A(i,j)=Pf;
       else 
           miula=kangla(i,j);
           Z5=miula/xigema33-1;
           Z6=(((1/xigema33)^2)*(xigemaxigema(i,j)^2))^0.5;
           baita3=Z5/Z6;
           Pf=1-normcdf(baita3,0,1);
           A(i,j)=Pf;
       end
    end
end
for i=322:521
    for j=1:2001
        xxx=0.25*(j-1);yyy=0.25*(i-1);
           xigema11=Kfind(diyi,xxx,yyy);xigema33=Kfind(disan3,xxx,yyy);%第二层含水层
         if xigema33>0
          miuc=c(i,j)*1e6;miufai=fai(i,j)/180;
          miuK=K((i-321),j);
          ZZZ=(miuK*miuc*cos(miufai)+sin(miufai)*(xigema11+xigema33)/2)/((xigema11-xigema33)/2)-1;
          ZZZZ=(((2*miuK*cos(miufai)/(xigema11-xigema33))^2)*(xigemac^2)+((xigema11+xigema33)/(((xigema11-xigema33))*cos(miufai)-2*miuK*miuc*sin(miufai)/(xigema11-xigema33))^2)*(xigemafai^2)+((2*miuc*cos(miufai)/(xigema11-xigema33))^2)*(xigemaK^2))^0.5;
          baita2=ZZZ/ZZZZ;
          Pf=1-normcdf(baita2,0,1);
          A(i,j)=Pf;
         else
             miula=kangla(i,j);
             Z7=miuK*miula/xigema33-1;
             Z8=(((miuK/xigema33)^2)*(xigemaxigema(i,j)^2)+((miula/xigema33)^2)*(xigemaK((i-321),j)^2))^0.5;
             baita4=Z7/Z8;
             Pf=1-normcdf(baita4,0,1);
             A(i,j)=Pf;
         end
    end
end
for i=522:721
    for j=1:2001            %第三层煤层
       xxx=0.25*(j-1);yyy=0.25*(i-1);
           xigema11=Kfind(diyi,xxx,yyy);xigema33=Kfind(disan3,xxx,yyy);
       if xigema33>0
        miuc=c(i,j)*1e6;miufai=fai(i,j)/180;
        Z=(miuc*cos(miufai)+sin(miufai)*(xigema11+xigema33)/2)/((xigema11-xigema33)/2)-1;
        ZZ=(((2*cos(miufai)/(xigema11-xigema33))^2)*(xigemac^2)+((xigema11+xigema33)/(((xigema11-xigema33))*cos(miufai)-2*miuc*sin(miufai)/(xigema11-xigema33))^2)*(xigemafai^2))^0.5;
        baita=Z/ZZ;
        Pf=1-normcdf(baita,0,1);
        A(i,j)=Pf;
       else 
           miula=kangla(i,j);
           Z5=miula/xigema33-1;
           Z6=(((1/xigema33)^2)*(xigemaxigema(i,j)^2))^0.5;
           baita3=Z5/Z6;
           Pf=1-normcdf(baita3,0,1);
           A(i,j)=Pf;
       end
    end
end
RE=A;
for p=649:681
    for q=601:1401
        RE(p,q)=1;
    end
end
