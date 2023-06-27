format long;
format compact;
clear all;
hold on;

f = '@(x) x';
p = '@(x) 1';
q = '@(x) -2';
nl0 = 1;
nl1 = -2;
nl = 0;
nr0 = 1;
nr1 = -1;
nr = 0;
l = -1;
r = 1;


solver(l,r,f,p,q,nl0,nl1,nl,nr0,nr1,nr);

%
function [X,Y] = solver(l,r,ft,pt,qt,nl0,nl1,nl,nr0,nr1,nr)
    
    %1 左边界,右边界,左边界初值,右边界初值,ui,节点个数
    %2 父节点
    %3 第一个子节点
    %4 区间左边界
    %5 区间右边界
    %6,7 α
    %8,9 β
    %10,11 δ
    %12,13,14 λ
    tree(1,1) = l;
    tree(1,2) = r;
    tree(1,3) = nl0;
    tree(1,4) = nl1;
    tree(1,5) = nl;
    tree(1,6) = nr0;
    tree(1,7) = nr1;
    tree(1,8) = nr;
    tree(1,11) = 1;
    tree(2,1) = 0;
    tree(3,1) = 0;
    tree(4,1) = l;
    tree(5,1) = r;
    tree(14,1) = 1;
    tree;

    %nl0(la+b)+nl1*a=nl
    %nr0(ra+b)+nr1*a=nr
    % a = (nl*nr0-nr*nl0)/(nl0*nr0(l-r)+nl1*nr0-nr1*nl0)
    % b = (nl-nl1*a)/nl0 - l*a
    tree(1,9) = (nl*nr0-nr*nl0)/(nl0*nr0*(l-r)+nl1*nr0-nr1*nl0);
    tree(1,10) = (nl-nl1*tree(1,9))/nl0 - l*tree(1,9);
   
    maxN = 10;


    %计算解并拓展二叉树
    flag = 0;
    C = 4;
    tree = upward(tree,ft,pt,qt);
    tree = downward(tree,ft,pt,qt);
    while flag == 0
        S = zeros(1,tree(1,11));
        clear sigma sigmaK sigma1 sigmaK1 old;
        len = 1;
        Sdiv = 0;
        [sigma,sigmaK] = check(tree,ft,pt,qt);
        for i = 1:tree(1,11)
            if tree(3,i) == 0
                S(i) = abs(sigmaK(i,maxN-1)) + abs(sigmaK(i,maxN)-sigmaK(i,maxN-2));
                if S(i) > Sdiv
                    Sdiv = S(i);
                end
            end
        end
        Sdiv = Sdiv/(2^C);
        for i = 1:tree(1,11)
            if S(i) > Sdiv
                l = tree(4,i);
                r = tree(5,i);
                lc = tree(1,11)+1;
                rc = lc+1;
                tree(3,i) = lc;
                tree(2,lc) = i;
                tree(2,rc) = i;
                tree(4,lc) = l;
                tree(5,rc) = r;
                tree(5,lc) = (l+r)/2;
                tree(4,rc) = (l+r)/2;
                tree(1,11) = tree(1,11)+2;
                old(len) = i;
                len = len + 1;
            end
        end
        tree = upward(tree,ft,pt,qt);
        tree = downward(tree,ft,pt,qt);
        [sigma1,sigmaK1] = check(tree,ft,pt,qt);
        error = 0;
        Max = 0;
        for it = 1:len-1
            i = old(it);
            %取2*maxM个等间距点，比较
            maxM = 10;
            s1 = getCheb(sigmaK(i,1:maxN+1),2*maxM);
            s2(1:maxM) = getCheb(sigmaK(i,1:maxN+1),maxM);
            s2(maxM+1:maxM*2) = getCheb(sigmaK(i,1:maxM+1),maxM);
            for j = 1:2*maxM
                error = max(error,abs(s1(j)-s2(j)));
                Max = max(Max,abs(s2(j)));
            end
        end
        if error/Max < 10^(-8) | tree(1,11)>100
            flag = 1;
        end
        %draw(tree,sigma1);
    end
    
    draw(tree,sigma1);
end

%绘制图像
function draw(tree,sigma)
    maxN = 10;
    persistent Ir;
    persistent Il;
    if isempty(Ir)
        [Il,Ir] = getInt(maxN);
    end
    for i = 1:tree(1,11)
        if tree(3,i) == 0
            nl = tree(4,i);
            nr = tree(5,i);
            for j = 1:maxN+1
                tau(j) = nl+(1+cos(2*pi*(j-1)/(2*maxN+1)))*(nr-nl)/2;
                [gl(j),gr(j),gl1(j),gr1(j),t] = g(tree,tau(j));
            end
            for j = 1:maxN+1
                s = gl(j)*gr1(j)-gl1(j)*gr(j);
                G(j,1:j) = gr(j)/s*gl(1:j);
                G(j,j:maxN+1) = gl(j)/s*gr(j:maxN+1);
            end
            for x = 1:maxN+1
                t = Il*(G(x,:).*sigma(i))';
                Y(x) = t(x);
            end
            plot(tau,Y,'b');
        end
    end
end

%计算每个节点的sigma
function [sigma,sigmaK] = check(tree,ft,pt,qt)
    maxN = 10;
    sigma = zeros(tree(1,11),2*maxN+1);
    sigmaK = zeros(tree(1,11),2*maxN+1);
    for i = 1:tree(1,11)
        %if tree(3,i) == 0
            nl = tree(4,i);
            nr = tree(5,i);
            P = discret(tree,i,pt,qt);
            for j = 1:maxN+1
                tau(j) = nl+(1+cos(2*pi*(j-1)/(2*maxN+1)))*(nr-nl)/2;
                [phil(j),phir(j)] = phi(tree,tau(j),pt,qt);
                tf(j) = getf(tree,tau(j),ft,pt,qt);
            end
            sigma(i,1:maxN+1) = pinv(P)*tf' + tree(12,i)*pinv(P)*phil' + tree(13,i)*pinv(P)*phir';
            for j = 1:maxN
               sigma(i,2*maxN+2-j) = sigma(i,j+1);
            end
            sigmaK(i,1:2*maxN+1) = ifft(sigma(i,1:2*maxN+1));
        %end
    end
end

%直接用P的离散化计算sigma，即1个节点的情况 （未使用）
function ret = direct(tree,ft,pt,qt)
    P = discret(tree,1,pt,qt)
    tau = l+(1+cos(linspace(0,2*pi,202)))*(r-l)/2;
    tau = tau(1:maxN+1);
    ft = zeros(1,maxN+1);
    for i = 1:maxN
        ft(i) = f(tau(i));
    end
    sigma = P\ft';
    for i = 1:maxN+1
        [gl(i),gr(i),gl1(i),gr1(i),t] = g(tree,tau(i));
    end
    for i = 1:maxN+1
        s = gl(i)*gr1(i)-gl1(i)*gr(i);
        G(i,1:i) = gr(i)/s*gl(1:i);
        G(i,i:maxN+1) = gl(i)/s*gr(i:maxN+1);
    end
    ret = G*sigma;
end


%计算α，β，δ
function ret = upward(tree,ft,pt,qt)
    maxN = 10;
    for it = 0:tree(1,11)-1
        i = tree(1,11)-it;
        if tree(3,i) == 0
            P = discret(tree,i,pt,qt);
            nl = tree(4,i);
            nr = tree(5,i);
            for j = 1:maxN+1
                
                tau(j) = nl+(1+cos(2*pi*(j-1)/(2*maxN+1)))*(nr-nl)/2;
                [phil(j),phir(j)] = phi(tree,tau(j),pt,qt);
                [gl(j),gr(j),t,t,t] = g(tree,tau(j));
                f(j) = getf(tree,tau(j),ft,pt,qt);
            end
            tree(6,i) = Inte(gl'.*(P\phil'))*(nr-nl)/2; %αl
            tree(7,i) = Inte(gr'.*(P\phil'))*(nr-nl)/2; %
            tree(8,i) = Inte(gl'.*(P\phir'))*(nr-nl)/2; %
            tree(9,i) = Inte(gr'.*(P\phir'))*(nr-nl)/2;
            tree(10,i) = Inte(gl'.*(P\f'))*(nr-nl)/2;
            tree(11,i) = Inte(gr'.*(P\f'))*(nr-nl)/2;
        else
            D = tree(3,i);
            E = tree(3,i)+1;
            alD = tree(6,D);
            alE = tree(6,E);
            arD = tree(7,D);
            arE = tree(7,E);
            blD = tree(8,D);
            blE = tree(8,E);
            brD = tree(9,D);
            brE = tree(9,E);
            dlD = tree(10,D);
            dlE = tree(10,E);
            drD = tree(11,D);
            drE = tree(11,E);
            delta = 1-arE*blD;
            tree(6,i) = (1-alE)*(alD-blD*arE)/delta + alE;
            tree(7,i) = arE*(1-brD)*(1-alD)/delta + arD;
            tree(8,i) = blD*(1-brE)*(1-alE)/delta + blE;
            tree(9,i) =(1-brD)*(brE-blD*arE)/delta + brD;
            tree(10,i) = (1-alE)*(dlD-blD*drE)/delta + dlE;
            tree(11,i) = (1-brD)*(drE-arE*dlD)/delta + drD;
        end
    end
    ret = tree;
end

%计算λ
function ret = downward(tree,ft,pt,qt)
    for i = 2:tree(1,11)
        if tree(3,i) ~= 0
            D = tree(3,i);
            E = tree(3,i)+1;
            alD = tree(6,D);
            alE = tree(6,E);
            arD = tree(7,D);
            arE = tree(7,E);
            blD = tree(8,D);
            blE = tree(8,E);
            brD = tree(9,D);
            brE = tree(9,E);
            dlD = tree(10,D);
            dlE = tree(10,E);
            drD = tree(11,D);
            drE = tree(11,E);
            delta = 1-arE*blD;
            mu = tree(14,i);
            mul = tree(12,i);
            mur = tree(13,i);
            tree(14,D) = mu;
            tree(14,E) = mu;
            tree(12,D) = mul;
            tree(13,E) = mur;
            a = mur*(1-brE)-mu*drE;
            b = mul*(1-alD)-mu*dlD;
            tree(12,E) = (-blD*a + b)/delta;
            tree(13,D) = (a - arE*b)/delta;
        end
    end
    ret = tree;
end

%离散化Pk
function [ret,m] = discret(tree,k,pt,qt)
    maxN = 10;
    persistent Ir;
    persistent Il;
    if isempty(Ir)
        [Il,Ir] = getInt(maxN);
    end
    nl = tree(4,k);
    nr = tree(5,k);
    tau = zeros(1,maxN+1);
    for i = 1:maxN+1
        tau(i) = nl+(1+cos(2*pi*(i-1)/(2*maxN+1)))*(nr-nl)/2;
        [phil(i),phir(i)] = phi(tree,tau(i),pt,qt);
        [gl(i),gr(i),t,t,t] = g(tree,tau(i));
    end
    ret = eye(maxN+1)+diag(phil)*Il*diag(gl)+diag(phir)*Ir*diag(gr);
    m = maxN+1;
end

%获取离散左、右积分的矩阵 （已检验）
function [Il,Ir] = getInt(n)
    x = zeros(1,2*n+1);
    Il = zeros(n+1,n+1);
    Ir = zeros(n+1,n+1);
    for i = 1:n+1
        x(i) = 1;
        if i>1
            x(2*n+3-i) = 1;
        end
        [t1,t2] = ChebInt(x);
        Il(:,i) = t1(1:n+1)';
        Ir(:,i) = t2(1:n+1)';
        x(i) = 0;
        if i>1
            x(2*n+3-i) = 0;
        end
    end
end

%获取Chebyshev函数的值
function ret = getCheb(x,n)
    m = length(x);
    tau = linspace(-1,1,n+1);
    ret = zeros(1,n);
    for i = 1:n
        for j = 1:m
            ret(i) = ret(i) + x(j)* cos((j-1)*acos(tau(i)));
        end
    end
end

%Chebyshev quadrature
function I = Inte(x)
    N = length(x);
    for i = 1:N-1
        x(2*N-i) = x(i+1);
    end
    a = ifft(x);
    
    I = 2*a(1);
    for i = 3:2*N-1
        if mod(i,2)==1
            I = I + 4 * a(i)/(1-(i-1)*(i-1));
        end
    end

end


%通过Chebshev node离散化左、右积分 （已检验）
function [retl,retr] = ChebInt(x)
    f = 2*ifft(x);
    n = floor(length(f)/2);
    f(n+1) = 0;
    
    for i = 2:n
        a(i) = (f(i-1)-f(i+1))/(i-1)/2;
        b(i) = (f(i+1)-f(i-1))/(i-1)/2;
        if mod(i,2) == 1
            a(1) = a(1) - a(i);
        else
            a(1) = a(1) + a(i);
        end
        b(1) = b(1) - b(i);
    end
    for i = 1:n-1
        a(i+1) = a(i+1)/2;
        a(length(f)+1-i) = a(i+1);
        b(i+1) = b(i+1)/2;
        b(length(f)+1-i) = b(i+1);
    end
    retl = fft(a);
    retr = fft(b);
end

%计算phi
%gl = nl0(x-l)-nl1
%gr = nr0(x-r)-nr1
%s = gl*gr'-gl;*gr
%phil = (p~gr'+q~gr)/s
%phir = (p~gl'+q~gl)/s
function [phil,phir] = phi(tree,x,pt,qt)
    p = str2func(pt);
    q = str2func(qt);
    [gl,gr,gl1,gr1,q0] = g(tree,x);
    s = gl*gr1-gl1*gr;
    phil = (p(x)*gr1+(q(x)-q0)*gr)/s;
    phir = (p(x)*gl1+(q(x)-q0)*gl)/s;
end

%计算gl，gr
function [gl,gr,gl1,gr1,q0] = g(tree,x)
    l = tree(1,1);
    r = tree(1,2);
    nl0 = tree(1,3);
    nl1 = tree(1,4);
    nr0 = tree(1,6);
    nr1 = tree(1,7);
    if abs(nl0)<abs(nl1) & abs(nr0)<abs(nr1)
        gl = nl1*cosh(x-l)-nl0*sinh(x-l);
        gr = nr1*cosh(x-r)-nr0*sinh(x-r);
        gl1 = nl1*sinh(x-l)-nl0*cosh(x-l);
        gr1 = nr1*sinh(x-r)-nr0*cosh(x-r);
        q0 = -1;
    else
        gl = nl0*(x-l)-nl1;
        gr = nr0*(x-r)-nr1;
        gl1 = nl0;
        gr1 = nr0;
        q0 = 0;
    end
end

%计算\tilda{f}
function ret = getf(tree,x,ft,pt,qt)
    f = str2func(ft);
    p = str2func(pt);
    q = str2func(qt);
    ret = f(x) - p(x)*tree(1,9) - q(x)*(tree(1,9)*x+tree(1,10));
end