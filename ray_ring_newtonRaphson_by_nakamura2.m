clear all
close all

%音速不均一な系の設定
grid_num = 1024;

%水領域
rwat = 100.e-3/2;%使用されていない[2018-06-10 追記：竹内]

%組織領域[2018-06-10 追記：竹内]
rtis = 70.e-3/2;
xtis = 110.e-3/2;
ytis = 110.e-3/2;

%脂肪領域[2018-06-10 追記：竹内]
rfat = 20.e-3/2;
xfat = 110.e-3/2 + 20.e-3;
yfat = 110.e-3/2;

%腫瘤領域[2018-06-10 追記：竹内]
rlsn = 20.e-3/2;
xlsn = 110.e-3/2 - 20.e-3;
ylsn = 110.e-3/2;

xmax = 110.e-3;
ymax = 110.e-3;
%下のような座標設定をすることで
% １：y　　　　　　　　２：　　----------→行
%　↑　　　　　　　　　　　　　|
%　 |　　　　　　　　　　　　　|
%　 |　　　　　　　　　＝　　　|　
%　 |　　　　　　　　 　 　　　↓
% 　-------------→x　　　　　　列
%このような関係が成り立ち，マトリクスをimagescしたときに１のように見え，直観的．[2018-06-10 追記：竹内]
x1 = linspace(0,xmax,grid_num);
xx = repmat(x1,grid_num,1);
y1 = linspace(ymax,0,grid_num);
yy = repmat(y1',1,grid_num);

R1 = (xx-xlsn).^2+(yy-ylsn).^2;
lsnmat = (R1<rlsn^2);%腫瘤領域をlogical型で表現[2018-06-10 追記：竹内]
R2 = (xx-xfat).^2+(yy-yfat).^2;
fatmat = (R2<rfat^2);%脂肪領域をlogical型で表現[2018-06-10 追記：竹内]
R3 = (xx-xtis).^2+(yy-ytis).^2;
w1 = (R3<rtis^2);
w02 = w1-fatmat;%組織が存在しうる領域から脂肪領域を取り除いた領域：2[2018-06-10 追記：竹内]
w01 = w1-lsnmat;%組織が存在しうる領域から腫瘤領域を取り除いた領域：１[2018-06-10 追記：竹内]
tismat = w01.*w02;%１でありかつ２である領域：真の組織領域[2018-06-10 追記：竹内]

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
lx = linspace(0+grid_size/2,L-grid_size/2,grid_num);%grid_sizeを保持するために始点・終点をずらした．結果の描写だけに用いるためずらしても計算結果に影響しない[2018-06-10 追記：竹内]
%3x3の平均値フィルターをかけスム‐シング
h = ones(3,3)*1/9;
n2 = filter2(h,n);
n(2:grid_num-1,2:grid_num-1) = n2(2:grid_num-1,2:grid_num-1);%パディング処理をしなかったためのエラーを含まないように代入の範囲を削減[2018-06-10 追記：竹内]

%送信の設定
ch = 256;
% theta_rg = linspace(0,2*pi,ch);%重複する箇所が出るのでは？[2018-06-10 追記：竹内]
%[2018-06-10-修正 竹内]
theta_rg = linspace(0, ((ch-1)/ch)*2*pi, ch);%センサ位置角
r_rg = 100.e-3/2;%リングトランスデューサ半径
x_rg = r_rg*cos(theta_rg)+L/2;
y_rg = r_rg*sin(theta_rg)+L/2;
ds = grid_size/8;

for tr_count = 1:1 %送信素子の選択
    tr = [x_rg(tr_count) y_rg(tr_count)];
    for re_count = 31:10:251
        re = [x_rg(re_count) y_rg(re_count)];
        theta0 = pi+theta_rg(tr_count);%送信素子と対向した位置角
        rpnum = 1;%音線作成回数
        r = tr;%記入する文字数を削減するため？
        loop = 1;%なんのループ？→おそらく角度修正にむけた内挿的処理ループ
        theta1 = theta0;%保存用？
        while(1)%なんのループ？
            dtheta = pi/180/(loop+1);%やりたいこと：最初は1ラジアンから刻み角を用意して，到達地点と受信素子との距離が十分近くなかったら走査範囲を細かくする．
            % 音線作成
            while(1)%なんのループ？
                % ここで2周目が開始しているニュアンス
                x(rpnum) = r(1); %#ok<SAGROW> 繰り返しごとにサイズ可変:送信素子のｘ座標　このループ内で変化している．
                y(rpnum) = r(2); %#ok<SAGROW> 繰り返しごとにサイズ可変:送信素子のｙ座標　デバッグ用か？
                if rpnum>2 %後退差分スキーム
                    %dx,dy : the change of x and y
                    dx = x(rpnum)-x(rpnum-1);
                    dy = y(rpnum)-y(rpnum-1);
                else
                    dx = ds*cos(theta1);
                    dy = ds*sin(theta1);
                end
                ix = round(x(rpnum)/grid_size+1);%ループごとに変化している．切り上げを行っている．
                jy = round(y(rpnum)/grid_size+1);%音線構築ループの各ステップにおける音線上の点を示すグリッド番号
%                 if ix>1||ix<grid_num+1
                    %nx,ny : the partial difference of n
                nx = (n(ix+1,jy)-n(ix-1,jy))/2/grid_size;%音線をストップするトリガーは他にあるはず．コメントアウト
%                 else
%                     nx = (n(ix+1,jy)-n(ix+1,jy))/2/grid_size;
%                 end
%                 if jy>1||jy<grid_num+1
                ny = (n(ix,jy+1)-n(ix,jy-1))/2/grid_size;
%                 else
%                     ny = (n(ix,jy+1)-n(ix,jy+1))/2/grid_size;
%                 end
                n_nearest = n(ix,jy);%どちらかが上書きされる．
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
            %音線作成終了
            %推定した音線の長さを計算する
            r = tr;
            rpnum = 1;
            rdis1 = re_distance;
            if rdis1<2.e-4
                figure;
                imagesc(lx,lx,n');hold on
                caxis([0.9 1.1]);set(gca,'Ydir','Normal')
                plot(tr(1),tr(2),'*');plot(re(1),re(2),'+')
                plot(x,y,'k')
                pause(2);%[2018-06-10 追記：竹内]
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