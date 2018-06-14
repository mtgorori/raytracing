clear all
close all

%音速不均一な系の設定
grid_num = 1024;

rwat = 100.e-3/2;

rtis = 70.e-3/2;
xtis = 110.e-3/2;
ytis = 110.e-3/2;

rfat = 20.e-3/2;
xfat = 110.e-3/2 + 20.e-3;
yfat = 110.e-3/2;

rlsn = 20.e-3/2;
xlsn = 110.e-3/2 - 20.e-3;
ylsn = 110.e-3/2;

xmax = 110.e-3;
ymax = 110.e-3;
x1 = linspace(0,xmax,grid_num);
xx = repmat(x1,grid_num,1);
y1 = linspace(ymax,0,grid_num);
yy = repmat(y1',1,grid_num);

R1 = (xx-xlsn).^2+(yy-ylsn).^2;
lsnmat = (R1<rlsn^2);
R2 = (xx-xfat).^2+(yy-yfat).^2;
fatmat = (R2<rfat^2);
R3 = (xx-xmax/2).^2+(yy-ymax/2).^2;
w1 = (R3<rtis^2);
w02 = w1-fatmat;
w01 = w1-lsnmat;
tismat = w01.*w02;

watmat = (R3>(rtis)^2);

us_wat = 1540;
us_lsn = 1550;
us_tis = 1540;
us_fat = 1530;

I = us_lsn.*lsnmat+us_fat.*fatmat+us_tis.*tismat+us_wat.*watmat;

n = us_wat./I;%n: distribution of n (refraction rate)

%スクリーンの設定
L = 110.e-3;
grid_size = L/grid_num;
lx = linspace(0+grid_size/2,L-grid_size/2,grid_num);
%3x3の平均値フィルターをかけスム‐シング
h = ones(3,3)*1/9;
n2 = filter2(h,n);
n(2:grid_num-1,2:grid_num-1) = n2(2:grid_num-1,2:grid_num-1);

%送信の設定
ch = 256;
theta_rg = linspace(0,2*pi*(ch-1)/ch,ch);
r_rg = 100.e-3/2;%リングトランスデューサ半径
x_rg = r_rg*cos(theta_rg)+L/2;
y_rg = r_rg*sin(theta_rg)+L/2;
ds = grid_size/8;

for tr_count = 1:1
    tr = [x_rg(tr_count) y_rg(tr_count)];
    figure;
    for re_count = 1:10:ch
        re = [x_rg(re_count) y_rg(re_count)];
        theta0 = pi+theta_rg(tr_count);
        rpnum = 1;
        r = tr;
        loop = 1;
        theta1 = theta0;
        %for loop = 1:10
        while(1)
            dtheta = pi/180/(loop+1);
            while(1)
                dtheta = pi/180/(loop+1);
                x(rpnum) = r(1);
                y(rpnum) = r(2);
                if rpnum>2
                    %dx,dy : the change of x and y
                    dx = x(rpnum)-x(rpnum-1);
                    dy = y(rpnum)-y(rpnum-1);
                else
                    dx = ds*cos(theta1);
                    dy = ds*sin(theta1);
                end
                ix = round(x(rpnum)/grid_size+1);
                jy = round(y(rpnum)/grid_size+1);
                if ix>1||ix<grid_num+1
                    %nx,ny : the partial difference of n
                    nx = (n(ix+1,jy)-n(ix-1,jy))/2/grid_size;
                else
                    nx = (n(ix+1,jy)-n(ix+1,jy))/2/grid_size;
                end
                if jy>1||jy<grid_num+1
                    ny = (n(ix,jy+1)-n(ix,jy-1))/2/grid_size;
                else
                    ny = (n(ix,jy+1)-n(ix,jy+1))/2/grid_size;
                end
                n_nearest = n(ix,jy);
                n_inter = n_nearest;
                detx = x(rpnum)/grid_size+1-ix;
                dety = y(rpnum)/grid_size+1-jy;
                if detx>=0
                    ix2 = round(x(rpnum)/grid_size+1+0.5);
                else
                    ix2 = round(x(rpnum)/grid_size+1-0.5);
                end
                if dety>=0
                    jy2 = round(y(rpnum)/grid_size+1+0.5);
                else
                    jy2 = round(y(rpnum)/grid_size+1-0.5);
                end
                if ix2>0&&ix>0&&ix2<grid_num+1&&ix<grid_num+1
                    if jy>0&&jy2>0&&jy<grid_num+1&&jy2<grid_num+1
                        lx1 = abs(x(rpnum)/grid_size+1-ix);lx2 = abs(x(rpnum)/grid_size+1-ix2);
                        ly1 = abs(y(rpnum)/grid_size+1-jy);ly2 = abs(y(rpnum)/grid_size+1-jy2);
                        n_inter = n(ix,jy)*lx2*ly2+n(ix2,jy)*lx1*ly2+n(ix,jy2)*lx2*ly1+n(ix2,jy2)*lx1*ly1;
                    end
                end
                DS = sqrt((dx+1/2/n_inter*(nx-(nx*dx/ds)*dx/ds)*grid_size^2)^2+(dy+1/2/n_inter*(ny-(ny*dy/ds)*dy/ds)*grid_size^2)^2);
                dsx = (dx+1/2/n_inter*(nx-(nx*dx/ds)*dx/ds)*grid_size^2)/DS*ds;
                dsy = (dy+1/2/n_inter*(ny-(ny*dy/ds)*dy/ds)*grid_size^2)/DS*ds;
                r(1) = x(rpnum)+dsx;
                r(2) = y(rpnum)+dsy;
                %境界の設定（計算の終了条件）
                R = sqrt((r(1)-xmax/2).^2+(r(2)-xmax/2).^2);
                if R>rwat
                    re_distance = sqrt((r(1)-re(1))^2+(r(2)-re(2))^2);
                    break
                end
                rpnum = rpnum+1;
            end    
            r = tr;
            rpnum = 1;
            rdis1 = re_distance;
            if rdis1<2.e-4
                imagesc(lx,lx,n');hold on
                caxis([0.9 1.1]);set(gca,'Ydir','Normal')
                plot(tr(1),tr(2),'*');plot(re(1),re(2),'+')
                plot(x,y,'k')
                plot(x_rg,y_rg,'k')
                pause(1);
%                 filename = ['Image',int2str(re_count)];
%                 saveas(figure(re_count),filename,'bmp')
                A = ['x_traced',int2str(re_count),'=x;'];
                B = ['y_traced',int2str(re_count),'=y;'];
                eval(A);
                eval(B);
                break
            end
            clear x y
            theta2 = theta1+dtheta;
            while(1)
                x(rpnum) = r(1);
                y(rpnum) = r(2);
                if rpnum>2
                    %dx,dy : the change of x and y
                    dx = x(rpnum)-x(rpnum-1);
                    dy = y(rpnum)-y(rpnum-1);
                else
                    dx = ds*cos(theta2);
                    dy = ds*sin(theta2);
                end
                ix = round(x(rpnum)/grid_size+1);
                jy = round(y(rpnum)/grid_size+1);
                if ix>1||ix<grid_num+1
                    %nx,ny : the partial difference of n
                    nx = (n(ix+1,jy)-n(ix-1,jy))/2/grid_size;
                else
                    nx = (n(ix+1,jy)-n(ix+1,jy))/2/grid_size;
                end
                if jy>1||jy<grid_num+1
                    ny = (n(ix,jy+1)-n(ix,jy-1))/2/grid_size;
                else
                    ny = (n(ix,jy+1)-n(ix,jy+1))/2/grid_size;
                end
                n_nearest = n(ix,jy);
                n_inter = n_nearest;
                detx = x(rpnum)/grid_size+1-ix;
                dety = y(rpnum)/grid_size+1-jy;
                if detx>=0
                    ix2 = round(x(rpnum)/grid_size+1+0.5);
                else
                    ix2 = round(x(rpnum)/grid_size+1-0.5);
                end
                if dety>=0
                    jy2 = round(y(rpnum)/grid_size+1+0.5);
                else
                    jy2 = round(y(rpnum)/grid_size+1-0.5);
                end
                if ix2>0&&ix>0&&ix2<grid_num+1&&ix<grid_num+1
                    if jy>0&&jy2>0&&jy<grid_num+1&&jy2<grid_num+1
                        lx1 = abs(x(rpnum)/grid_size+1-ix);lx2 = abs(x(rpnum)/grid_size+1-ix2);
                        ly1 = abs(y(rpnum)/grid_size+1-jy);ly2 = abs(y(rpnum)/grid_size+1-jy2);
                        n_inter = n(ix,jy)*lx2*ly2+n(ix2,jy)*lx1*ly2+n(ix,jy2)*lx2*ly1+n(ix2,jy2)*lx1*ly1;
                    end
                end
                DS = sqrt((dx+1/2/n_inter*(nx-(nx*dx/ds)*dx/ds)*grid_size^2)^2+(dy+1/2/n_inter*(ny-(ny*dy/ds)*dy/ds)*grid_size^2)^2);
                dsx = (dx+1/2/n_inter*(nx-(nx*dx/ds)*dx/ds)*grid_size^2)/DS*ds;
                dsy = (dy+1/2/n_inter*(ny-(ny*dy/ds)*dy/ds)*grid_size^2)/DS*ds;
                r(1) = x(rpnum)+dsx;
                r(2) = y(rpnum)+dsy;
                %境界の設定（計算の終了条件）
                R = sqrt((r(1)-xmax/2).^2+(r(2)-xmax/2).^2);
                if R>rwat
                    re_distance = sqrt((r(1)-re(1))^2+(r(2)-re(2))^2);
                    break
                end
                rpnum = rpnum+1;
            end   
            rdis2 = re_distance;
            theta1 = theta1-rdis1*dtheta/(rdis2-rdis1);
            r = tr;
            rpnum = 1;
            z(loop) = theta1;
            u(loop) = rdis1;
            loop = loop+1;
            clear x y
        end
    end
end

figure;imagesc(lx,lx,n');hold on
caxis([0.9 1.1]);set(gca,'Ydir','Normal')
plot(tr(1),tr(2),'*')
plot(re(1),re(2),'+')
for re_count = 1:ch
    C = ['x = x_traced',int2str(re_count),';'];
    D = ['y = y_traced',int2str(re_count),';'];
    eval(C);
    eval(D);
    plot(x,y,'r')
end