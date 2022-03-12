// Définition des fonctions

// fct critère
function y=f(x)
    y=2*x
endfunction 

// gradient de la fct critère
function y=gradf(x)
    y=2
endfunction

// fct contraite
function y=g(x)
    y=1-x
endfunction 

 // gradient de la fct contrainte
function y=gradg(x)
    y=-1
endfunction

// Définition de la fct pénalisée avec barrière inverse
function z=p(x,r)
    if (g(x)) < 0 then 
        z=(f(x)-r*((g(x)).^(-1)));
    else 
        z=f(x);
    end
endfunction 

// grad de la fct a minimiser
function z=gradp(x,r)
    if (g(x)) < 0 then
        z=(gradf(x)+r*gradg(x)*((g(x)).^(-2))); 
    else 
        z=gradf(x);
    end
endfunction 

// parametre pour armijo
alphainit = 1;
beta = 0.001;
tau = 0.3;

// parametre pour descente gradient
eps = 10^(-4);
Niter = 100;

// Exo 1 Question b ----- Méthode du gradient pour la fonction pénalisée

function x=gradient(x,eps,Niter,r)
    x=x0
    k=0
    while (norm(gradp(x,r))>eps)&(k<Niter)&(g(x)<0)
        d=-gradp(x,r)
        alpha=alphainit;
        while p(x+alpha*d,r)>p(x,r)+beta*alpha*(d'*gradp(x,r)) 
            alpha=tau*alpha
        end
        x=x+alpha*d
        k=k+1
    end
endfunction

// Pour x0=2
x0=2; r=3; x=gradient(x0,eps,Niter,r);

// Exo 1 Question c
resid=norm(f(x0)-f(x));
disp(resid,"Résidu");

// Exo 1 Question d
mu=0.5; // taux de pénalisation de r

function [Nit,r,Xnew]=penaliteinter(x0,r,Niter,eps)
    xnew=x0;
    xold=0;
    k=0;
    Nit=0;
    while ((norm(f(xnew)-f(xold))>eps*norm(f(xold)))&(k<Niter)&(g(xnew)<0))
        xold=xnew;
        xnew=gradient(xnew,eps,Niter,r); // On utilise la méthode du gradient en initialisant au dernier point trouvé
        k=k+1;
        Nit=k;
        r=mu*r;
    end
    xnew;
    r;
    disp('valeur obtenue :');disp(xnew);
endfunction
[r,xnew]=penaliteinter(x0,r,Niter,eps);
disp(xnew)


/* ******** Decente grad projeté Ex 2 ***********

function [Nit,Xnew] = gradientPrejete(X0,Niter,eps)
    Xnew=X0;
    Xold = 0;
    k = 0;
    Nit = 0;
    while ((norm(f(Xnew) - f(Xold)) > eps * f(Xold)) & (k<Niter))
        Xold = Xnew;
        Xnew1=gradientf(Xnew,eps,Niter);
        Xnew=project(Xnew1)
        k = k+1;
        Nit = k;
    end
    disp('valeur obtenue');disp(Xnew)
endfunction
X0 = 2;
function P=project(X,a,b)
    [n,m]=size(X);
    for i=i:n
        P(i)=min(max(a(i),X(i)),b(i));
    end
    P;
endfunction

a = [-1;-%inf]; b = [+%inf;5]; X0=[-11;5];
P=project(X,a,b)

*/
