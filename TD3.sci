//Fichier qui contient toute les fonctions du TD3, fct de calcul de pas et fct de descente

//Fonctions qui interviennent dans le calcul de certain pas

// Utilisée pour le pas d'armijo
// Fonction f exo 4 // à refaire en fct de l'exercice
function res=f(x)
    x1=x(1,1)
    x2=x(2,1)
    res=5*x1^2+0.5*x2^2-3*(x1+x2)
endfunction

// Utilisée pour le pas d'armijo
// Calcul gradient fonction quadratique // à recoder (calcul du gradient à la main)
function [grad]=calculGradient(A,b,x)
    grad=-(A*x)-b
endfunction

// --------------------------------------------------------------------------------------------------------//
// Fonction des pas

// Calcul du pas optimal dans le cas d'une fonction quadratique 
function rho=pasCauchyQuadra(A,b,x,d)
    rho=-((A*x+b)'*d)/(d'*A*d);
endfunction

// Cas d'une fonction de moindres carrés
function rho=pasCauchyMc(A,b,x,d)
    rho=-((A'*(A*x-b))'*d)/((norm(A*d))^2)
endfunction

// Pas d'Armijo 
function rho=pasArmijo(x,d,omega,tau,rho0,A,b)
    rho=rho0
    r=calculGradient(A,b,x)
   // r=gradf(x) //calculer à la main 
    while (f(x+rho*d)>f(x)+omega*rho*r'*d)
        rho=tau*rho
    end
endfunction

// --------------------------------------------------------------------------------------------------------//

// Fonction de descente avec différents pas

// Exo 4 question 2
x0 = [-2 ; -7]
A = [10, 0; 0, 1]
b = [-3; -3]
nMax = 4
function [dk,rho,xk]=exoQuadra(A,b,x0,nMax)
    k=0
    xk=x0
    while ((k<nMax) | (k==0))
        disp("Itération k : ")
        disp(k + 1)
        dk=-(A*xk)-b
        disp("Direction de descente d_k-1: ")
        disp(dk)
        rho=pasCauchyQuadra(A,b,xk,dk)
        disp("Pas optimal rho : ")
        disp(rho)
        xk=xk+rho*dk
        disp("Itéré courant x_k: ")
        disp(xk)
        k=k+1
    end
endfunction

[D,B,C] = exoQuadra(A,b,x0,nMax)


// Exo 4 question 2
x0 = [-2 ; -7]
d0 = [23; 10]
omega = 10^-4
tau = 0.7
rho0 = 1
A = [10, 0; 0, 1]
b = [-3; -3]
function [rho,fonc,x1]=exoArmijo(x0,d0,omega,tau,rho0,A,b)
    //Calcul pas d'armijo 
    rho=rho0
    r=calculGradient(A,b,x0) 

    //Toute premiere étape i=0
    i=0
    while (f(x0+rho*d0)>f(x0)+omega*rho*r'*d0)
        
        disp("Etape i :")
        disp(i)
        if (i>0)
            rho=tau*rho
        end
        disp("Pas optimal rho : ")
        disp(rho)   
        //Affichage du gros calcul demandé exo 4 question 3 
        calcul=f(x0+rho*d0)-f(x0)-omega*rho*r'*d0
        disp("grosse opération sur f")
        disp(calcul)
        i=i+1
    end
    //Affichage de X1
    x1=x0+rho*d0
    disp("Itéré X1 : ")
    disp(x1)  
endfunction

//[D,B,C] = exoArmijo(x0,d0,omega,tau,rho0,A,b)
