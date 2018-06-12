clear all
close all

%�����s�ψ�Ȍn�̐ݒ�
grid_num = 1024;

%���̈�
rwat = 100.e-3/2;%�g�p����Ă��Ȃ�[2018-06-10 �ǋL�F�|��]

%�g�D�̈�[2018-06-10 �ǋL�F�|��]
rtis = 70.e-3/2;
xtis = 110.e-3/2;
ytis = 110.e-3/2;

%���b�̈�[2018-06-10 �ǋL�F�|��]
rfat = 20.e-3/2;
xfat = 110.e-3/2 + 20.e-3;
yfat = 110.e-3/2;

%��᎗̈�[2018-06-10 �ǋL�F�|��]
rlsn = 20.e-3/2;
xlsn = 110.e-3/2 - 20.e-3;
ylsn = 110.e-3/2;

xmax = 110.e-3;
ymax = 110.e-3;
%���̂悤�ȍ��W�ݒ�����邱�Ƃ�
% �P�Fy�@�@�@�@�@�@�@�@�Q�F�@�@----------���s
%�@���@�@�@�@�@�@�@�@�@�@�@�@�@|
%�@ |�@�@�@�@�@�@�@�@�@�@�@�@�@|
%�@ |�@�@�@�@�@�@�@�@�@���@�@�@|�@
%�@ |�@�@�@�@�@�@�@�@ �@ �@�@�@��
% �@-------------��x�@�@�@�@�@�@��
%���̂悤�Ȋ֌W�����藧���C�}�g���N�X��imagesc�����Ƃ��ɂP�̂悤�Ɍ����C���ϓI�D[2018-06-10 �ǋL�F�|��]
x1 = linspace(0,xmax,grid_num);
xx = repmat(x1,grid_num,1);
y1 = linspace(ymax,0,grid_num);
yy = repmat(y1',1,grid_num);

R1 = (xx-xlsn).^2+(yy-ylsn).^2;
lsnmat = (R1<rlsn^2);%��᎗̈��logical�^�ŕ\��[2018-06-10 �ǋL�F�|��]
R2 = (xx-xfat).^2+(yy-yfat).^2;
fatmat = (R2<rfat^2);%���b�̈��logical�^�ŕ\��[2018-06-10 �ǋL�F�|��]
R3 = (xx-xtis).^2+(yy-ytis).^2;
w1 = (R3<rtis^2);
w02 = w1-fatmat;%�g�D�����݂�����̈悩�玉�b�̈����菜�����̈�F2[2018-06-10 �ǋL�F�|��]
w01 = w1-lsnmat;%�g�D�����݂�����̈悩���᎗̈����菜�����̈�F�P[2018-06-10 �ǋL�F�|��]
tismat = w01.*w02;%�P�ł��肩�Q�ł���̈�F�^�̑g�D�̈�[2018-06-10 �ǋL�F�|��]

watmat = (R3>(rtis)^2);

us_wat = 1540;
us_lsn = 1550;
us_tis = 1540;
us_fat = 1530;

I = us_lsn.*lsnmat+us_fat.*fatmat+us_tis.*tismat+us_wat.*watmat;

n = us_wat./I;%n: distribution of n (refraction rate)

%�X�N���[���̐ݒ�
L = 110.e-3;
grid_size = L/grid_num;
lx = linspace(0+grid_size/2,L-grid_size/2,grid_num);%grid_size��ێ����邽�߂Ɏn�_�E�I�_�����炵���D���ʂ̕`�ʂ����ɗp���邽�߂��炵�Ă��v�Z���ʂɉe�����Ȃ�[2018-06-10 �ǋL�F�|��]
%3x3�̕��ϒl�t�B���^�[�������X���]�V���O
h = ones(3,3)*1/9;
n2 = filter2(h,n);
n(2:grid_num-1,2:grid_num-1) = n2(2:grid_num-1,2:grid_num-1);%�p�f�B���O���������Ȃ��������߂̃G���[���܂܂Ȃ��悤�ɑ���͈̔͂��팸[2018-06-10 �ǋL�F�|��]

%���M�̐ݒ�
ch = 256;
% theta_rg = linspace(0,2*pi,ch);%�d������ӏ����o��̂ł́H[2018-06-10 �ǋL�F�|��]
%[2018-06-10-�C�� �|��]
theta_rg = linspace(0, ((ch-1)/ch)*2*pi, ch);%�Z���T�ʒu�p
r_rg = 100.e-3/2;%�����O�g�����X�f���[�T���a
x_rg = r_rg*cos(theta_rg)+L/2;
y_rg = r_rg*sin(theta_rg)+L/2;
ds = grid_size/8;

for tr_count = 1:1 %���M�f�q�̑I��
    tr = [x_rg(tr_count) y_rg(tr_count)];
    for re_count = 31:10:251
        re = [x_rg(re_count) y_rg(re_count)];
        theta0 = pi+theta_rg(tr_count);%���M�f�q�ƑΌ������ʒu�p
        rpnum = 1;%�����쐬��
        r = tr;%�L�����镶�������팸���邽�߁H
        loop = 1;%�Ȃ�̃��[�v�H�������炭�p�x�C���ɂނ������}�I�������[�v
        theta1 = theta0;%�ۑ��p�H
        while(1)%�Ȃ�̃��[�v�H
            dtheta = pi/180/(loop+1);%��肽�����ƁF�ŏ���1���W�A�����獏�݊p��p�ӂ��āC���B�n�_�Ǝ�M�f�q�Ƃ̋������\���߂��Ȃ������瑖���͈͂��ׂ�������D
            % �����쐬
            while(1)%�Ȃ�̃��[�v�H
                % ������2���ڂ��J�n���Ă���j���A���X
                x(rpnum) = r(1); %#ok<SAGROW> �J��Ԃ����ƂɃT�C�Y��:���M�f�q�̂����W�@���̃��[�v���ŕω����Ă���D
                y(rpnum) = r(2); %#ok<SAGROW> �J��Ԃ����ƂɃT�C�Y��:���M�f�q�̂����W�@�f�o�b�O�p���H
                if rpnum>2 %��ލ����X�L�[��
                    %dx,dy : the change of x and y
                    dx = x(rpnum)-x(rpnum-1);
                    dy = y(rpnum)-y(rpnum-1);
                else
                    dx = ds*cos(theta1);
                    dy = ds*sin(theta1);
                end
                ix = round(x(rpnum)/grid_size+1);%���[�v���Ƃɕω����Ă���D�؂�グ���s���Ă���D
                jy = round(y(rpnum)/grid_size+1);%�����\�z���[�v�̊e�X�e�b�v�ɂ����鉹����̓_�������O���b�h�ԍ�
%                 if ix>1||ix<grid_num+1
                    %nx,ny : the partial difference of n
                nx = (n(ix+1,jy)-n(ix-1,jy))/2/grid_size;%�������X�g�b�v����g���K�[�͑��ɂ���͂��D�R�����g�A�E�g
%                 else
%                     nx = (n(ix+1,jy)-n(ix+1,jy))/2/grid_size;
%                 end
%                 if jy>1||jy<grid_num+1
                ny = (n(ix,jy+1)-n(ix,jy-1))/2/grid_size;
%                 else
%                     ny = (n(ix,jy+1)-n(ix,jy+1))/2/grid_size;
%                 end
                n_nearest = n(ix,jy);%�ǂ��炩���㏑�������D
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
                %���E�̐ݒ�i�v�Z�̏I�������j
                R = sqrt((r(1)-xmax/2).^2+(r(2)-xmax/2).^2);
                if R>rwat
                    re_distance = sqrt((r(1)-re(1))^2+(r(2)-re(2))^2);
                    break
                end
                rpnum = rpnum+1;
            end
            %�����쐬�I��
            %���肵�������̒������v�Z����
            r = tr;
            rpnum = 1;
            rdis1 = re_distance;
            if rdis1<2.e-4
                figure;
                imagesc(lx,lx,n');hold on
                caxis([0.9 1.1]);set(gca,'Ydir','Normal')
                plot(tr(1),tr(2),'*');plot(re(1),re(2),'+')
                plot(x,y,'k')
                pause(2);%[2018-06-10 �ǋL�F�|��]
                close;
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
                %���E�̐ݒ�i�v�Z�̏I�������j
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