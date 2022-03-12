// FONCTION QUI AIDE POUR LES ALGOS EN DESSOUS

// f est une fct quadratique sous sa forme général
function [y]=f(A,b,c,x)
    y= c + b'*x+(1/2)*x'*A*x
endfunction

// fonction du calcul du gradient pour une fct quadratique
function r=gradf(x,A,b)
    r=b+A*x
endfunction
//----------------------------------------------------------//
// f est une fct à implémenter a chaque fois
function [y]=fNQ(x)
    x1=x(1,1)
    x2=x(2,1)
    y=x1^4+4*x2^4+4*x1*x2
endfunction

// fonction du calcul du gradient à implémenter à chaque fois
function r=gradfNQ(x)
    x1=x(1,1)
    x2=x(2,1)
    r=[4*x1^3+4*x2 ; 16*x2^3+4*x1]
endfunction

// fonction du calcul du hessien à implémenter à chaque fois
// hessien doit être inversible
function H=hessienNQ(x)
    x1=x(1,1)
    x2=x(2,1)
    H=[12*x1^2 4;4 42*x2^2]
endfunction

// fonction du produit scalaire
function s=produitScalaire(u,v)
    s=0.5*(norm(u)^2+norm(v)^2-norm(v-u)^2);
endfunction

//------------------------------------------------------------------------------//
//IMPLÉMENTATION DE L'ALGO : GRADIENT CONJUGUÉ AVEC PAS OPTIMAL DE CAUCHY POUR FCT QUADRATIQUE

// Calcul du pas optimal de cauchy dans le cas d'une fonction quadratique 
function rho=pasCauchyQuadra(A,b,x,d)
    rho=-((A*x+b)'*d)/(d'*A*d);
endfunction

// calcul de Beta pour gradient conjugué avec fct quadratique
function B=calculeBeta(x,x1,A,b)
    B=(norm(gradf(x1,A,b))*norm(gradf(x1,A,b)))/(norm(gradf(x,A,b)*norm(gradf(x,A,b))))
end

// Algorithme du gradient conjugué pour les fcts quadratiques
function x=gradientConjuguePasCauchyFctQuadra(x0,epsi,A,b,c)
    k=0;
    x=x0;
    d=-A*x0-b;
    while (norm(gradf(x,A,b)) > epsi) 
        //rho= (norm(gradf(x,A,b))*norm(gradf(x,A,b))) / (norm(A*d)*norm(d));
        rho=pasCauchyQuadra(A,b,x,d);
        x1=x+rho*d;
        B1=calculeBeta(x,x1,A,b);
        d=-gradf(x1,A,b)+B1*d;
        x=x1;
        k=k+1;
    end
    disp("nb itération")
    disp(k)
endfunction

// tester le gradient conjugué pour la fct de l'exo3 td4


//------------------------------------------------------------------------------//
// IMPLÉMENTATION DE L'ALGO : GRADIENT CONJUGUÉ AVEC PAS ARMIJO VARIANTE FLECHTER-REEVES/POLAK RIBIÈRE AVEC RELANCE
//

// Pas d'Armijo 
function rho=pasArmijo(x,d,omega,tau,rho0)
    rho=rho0
    r=gradfNQ(x)
   // r=gradf(x) //calculer à la main 
    while (fNQ(x+rho*d)>fNQ(x)+omega*rho*r'*d)
        rho=tau*rho
    end
endfunction

// Flechter-Reeves
function B=flechterReeves(x1,x)
    numerateur=norm(gradfNQ(x1))^2;
    denominateur=norm(gradfNQ(x))^2;
    B=numerateur/denominateur;
endfunction

// Polak Ribière
function B=polakRibiere(x1,x)
    num1=gradfNQ(x1);
    num2=gradfNQ(x1)-gradfNQ(x);
    numerateur=produitScalaire(num1,num2);
    denominateur=gradfNQ(x)^2;
    B=numerateur/denominateur;
endfunction

// Algo gradient conjugué fct non quadratique avec relance
function x=gradientConjugueFctNQRelance(x0,epsi,nbrelance)
    k=0;
    x=x0;
    d=-gradfNQ(x0);
    rho=1;
    while (norm(gradfNQ(x)) > epsi)
        rho=pasArmijo(x,d,10^(-4),10^(-2),rho);
        x1=x+rho*d;
        // On choisit entre Flechter-Reeves et Polak Ribière
        B=flechterReeves(x1,x);
        d=-gradfNQ(x1)+B*d;
        k=k+1;
        x=x1;
        if (k == nbrelance) then
            d=-gradfNQ(x);
            k=0;
        end
    end
endfunction


//------------------------------------------------------------------------------//
// IMPLÉMENTATION DE L'ALGO : DES ALGOS DE NEWTON ET BFGS

// Algo de Newton
function x=newtonNQ(x0,epsi,nbmax)
    k=0
    x=x0
    d=-inv(hessienNQ(x))*gradfNQ(x)
    disp(produitScalaire(gradfNQ(x),inv(hessienNQ(x))*gradfNQ(x)))
    // problème sur le while condition supérieur largement à epsilon
    while ((produitScalaire(gradfNQ(x),inv(hessienNQ(x))*gradfNQ(x)) < epsi) || (k > nbmax))
        x=x+d
        k=k+1
    end
disp(k)
endfunction

// fonction de calcul de Sk+1 par rapport a Sk
function Sk1=approxHessien(Sk,y,delta)
    num1=y'*Sk*y
    den1=delta'*y
    frac1=num1/den1
    num2=delta*delta'
    den2=delta'*y
    frac2=num2/den2
    num3=delta*y'*Sk+Sk*y*delta'
    den3=delta'*y
    frac3=num3/den3
    Sk1=Sk+(1+frac1)*frac2-frac3
endfunction

// Algo de BFGS
// On prend S0 symétrique et définie positive
function x=BFGSNQ(x0,epsi)
    k=0
    x=x0
    S=[1 0; 0 1]
    rho=1
    disp(produitScalaire(gradfNQ(x),inv(hessienNQ(x))*gradfNQ(x)))
    while ((produitScalaire(gradfNQ(x),inv(hessienNQ(x))*gradfNQ(x))) < epsi)
        d=-S*gradfNQ(x)
        rho=pasArmijo(x,d,10^(-4),10^(-2),rho)
        x1=x+rho*d
        delta=x1-x
        y=gradfNQ(x1)-gradfNQ(x)
        S=approxHessien(S,y,delta)
        x=x1
        k=k+1
    end
disp(k)
endfunction
