# -*- coding: utf-8 -*-
"""
Created on Wed May 12 15:33:02 2021

@author: elise.foulatier
"""

import random
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

### DANS CE FICHIER ON CHERCHE A SIMULER LA PROPAGATION DU FEU EN DONNANT EN ENTREE
### LES VITESSES DU VENT SUIVANT X ET Y

### PARAMETRES D'ENTREE
#vitesse suivant x et y en km/h
Vx = 60
Vy = 60

#on déduit la norme du vent
V=sqrt(Vx**2+Vy**2)

#Vitesse à partir de laquelle le feu ne peut pas se propager dans le sens opposé à celui du vent
Vmax=60


###VALEURS DES PROBAS P1 ET P2 REGISSANT LA SIMULATION DE LA PROPAGATION DU FEU

#rappel des coefficients du modèle polynomial obtenu précédemment
#coefficients du polynome du second degré
alpha0 = 4.52790762e-01
alpha1 = 9.58264437e-04
alpha2 = 3.61499382e-05

#proba p1
if V<Vmax:
    p1 = alpha0 + alpha1*V + alpha2*(V**2)
else :
    p1 = alpha0 + alpha1*Vmax + alpha2*(Vmax**2)
    
p2 = 0.3

##################################################################################
### DEFINITION D'UN CERTAIN NOMBRE DE FONCTIONS UTILES AU COURS DE LA SIMULATION
##################################################################################

#fonction prenant en paramètre une matrice 
#retournant toutes les coordonnées des points rouges (coefficients égaux à 1)
def casesRouges(Matrice):
    M=len(Matrice[0]) #nombre d'abscisses
    N=len(Matrice) #nombre d'ordonnées (a priori la matrice n'est pas forcément carrée)
    
    ListePointsRouges=[] #liste qui contiendra les coordonées de tous les points rouges
    
    #balayage de la matrice
    for i in range(N):
        for j in range(M):
            if Matrice[i][j]==1:
                ListePointsRouges.append([j,i])  #on ajoute les coordonnées du point rouge rencontré à ListePointsRouges
    return ListePointsRouges


#fonction qui permet l'affichage de l'évolution du feu 
#paramètre d'entrée : matrice modélisant l'espace 
#retourne 6 listes: 
    #les points verts (leurs abscisses et leurs ordonnées) ;
    #les points rouges (leurs abscisses et leurs ordonnées) ; 
    #les points noirs (leurs abscisses et leurs ordonnées)
def affichageCases(Matrice):
    #nombre d'abscisses
    M=len(Matrice[0])
    #nombre d'ordonnées
    N=len(Matrice)
    Xverts=[]
    Yverts=[]
    Xrouges=[]
    Yrouges=[]
    Xnoirs=[]
    Ynoirs=[]
    #balayage de la matrice
    for i in range(N):
        for k in range(M):
            #cas des points verts
            if Matrice[i][k]==0:
                Xverts.append(longueur*k) #on retient l'abscisse correspondante
                Yverts.append(longueur*i) #on retient l'ordonnee correspondante
            #cas des points rouges
            elif Matrice[i][k]==1:
                Xrouges.append(longueur*k) #on retient l'abscisse correspondante
                Yrouges.append(longueur*i) #on retient l'ordonnee correspondante
            #cas des points noirs
            else :
                Xnoirs.append(longueur*k) #on retient l'abscisse correspondante
                Ynoirs.append(longueur*i) #on retient l'ordonnee correspondante
            
    return Xverts,Yverts,Xrouges,Yrouges,Xnoirs,Ynoirs

### DEFINITION DE L'ESPACE
#Longueur totale d'un côté de l'espace carré modélisé (en m)
L = 1000
#Espace entre deux points d'observation (en m)
longueur = 50 
#nombre de points observés sur une ligne du quadrillage
N = int(L/longueur) + 1 
#représentation du quadrillage par une matrice carrée de taille N*N (modèle bidimensionnel)
Matrice = np.zeros((N,N))

#AFFICHAGE DU MAILLAGE AVANT DEPART DU FEU (que des points verts)
#Construction de la liste des abscisses
X=[]
for i in range(N):
    X.append(longueur*i)
    
#Affichage du maillage ligne par ligne avant départ du feu
for i in range(N):
    plt.plot(X,[i*longueur]*N,"o",color="green")
    
#Position initiale du feu (les coordonnées du point doivent être multiples de longueur = 50)
#coordonnées en m (position réelle du point de départ du feu)
x0 = 500
y0= 500
#Couple d'entiers désignant la position du point dans Matrice
m0=int(x0/longueur)
n0=int(y0/longueur)
Matrice[n0][m0]=1
plt.plot([x0],[y0],"o",color="red")
#Affichage de la flèche désignant l'orientation du vent
plt.quiver(L+longueur,L+longueur,Vx,Vy,color="blue",label="direction du vent")
plt.legend()
#pause d'une seconde pour avoir le temps d'identifier sur la figure la position du point de départ
plt.pause(1)

#DEFINITION DU POINT D'ARRET DU FEU
xfinal=x0
yfinal=y0
    
#tant qu'aucune extrémité n'est atteinte on augmente la position du point d'arrêt suivant x et y (en fonction de Vx et Vy)
while not (xfinal>L or xfinal<0 or yfinal>L or yfinal<0):
    xfinal += longueur*Vx*0.01
    yfinal += longueur*Vy*0.01
    
#on obtient ainsi le point pour lequel on arrête la simulation une fois qu'il est atteint
if xfinal>L:
    entierMax = [N-1,int(yfinal//longueur)]
elif xfinal<0 :
    entierMax = [0,int(yfinal//longueur)]
elif yfinal>L:
    entierMax = [int(xfinal//longueur),N-1]
else :
    entierMax = [int(xfinal//longueur),0]


### Affectation de coefficients multiplicateurs pour modifier la probabilité en fonction de la direction de propagation

if Vx>0:
    alphaEstOuest=abs(Vx/V)+1
    alphaOuestEst=1-abs(Vx/V)
else :
    alphaOuestEst=abs(Vx/V)+1
    alphaEstOuest=1-abs(Vx/V)
    
if Vy>0:
    alphaSudNord=abs(Vy/V)+1
    alphaNordSud=1-abs(Vy/V)
else:
    alphaNordSud=abs(Vy/V)+1
    alphaSudNord=1-abs(Vy/V)
    

j=False
k=0
while casesRouges(Matrice)!=[]:
        ListeCasesRouges = casesRouges(Matrice)
        #on parcourt la liste des points rouges et on envisage pour chacun d'eux la propagation du feu dans 4 directions
        for caseRouge in ListeCasesRouges :
            #coordonnées du point brulant considéré
            #entier désignant l'abscisse
            m=caseRouge[0]
            #entier désignant l'ordonnée
            n=caseRouge[1]
        
            if [m,n]==entierMax: #si un bord du quadrillage est atteint j=Vrai et la simulation de la propagation s'arrête
                j=True
        
            if m<N-1 and j==False :
                if Matrice[n][m+1]==0:
                    #propagation du feu dans la direction x croissants
                    nombre_aleat_case_voisine=random.randint(1,100)
                    if nombre_aleat_case_voisine<alphaEstOuest*p1*100 :
                        Matrice[n][m+1]=1
        
            if m>0 and j==False:
                if Matrice[n][m-1]==0:
                    #propagation du feu dans la direction x décroissants
                    nombre_aleat_case_voisine=random.randint(1,100)
                    if nombre_aleat_case_voisine<alphaOuestEst*p1*100 :
                        Matrice[n][m-1]=1
                
            if n<N-1 and j==False:
                if Matrice[n+1][m]==0:
                    #propagation du feu dans la direction y croissants
                    nombre_aleat_case_voisine=random.randint(1,100)
                    if nombre_aleat_case_voisine<alphaSudNord*p1*100 :
                        Matrice[n+1][m]=1
        
            if n>0 and j==False:
                if Matrice[n-1][m]==0:
                    #propagation du feu dans la direction y décroissants
                    nombre_aleat_case_voisine=random.randint(1,100)
                    if nombre_aleat_case_voisine<alphaNordSud*p1*100 : 
                        Matrice[n-1][m]=1         
            #évolution de l'état d'une case rouge
            nombre_aleat_reste_brulant = random.randint(1,100)
            if nombre_aleat_reste_brulant < p2*100: #si le nombre généré est inférieur à 100 fois la probabilité alors la case passe à l'état brûlé
                Matrice[n][m]=2                     #la case considérée est désormais brûlée, coef égal à 2
            
        
        #AFFICHAGE DES POINTS DANS MATPLOTLIB
        Xverts=affichageCases(Matrice)[0]
        Yverts=affichageCases(Matrice)[1]
        Xrouges=affichageCases(Matrice)[2]
        Yrouges=affichageCases(Matrice)[3]
        Xnoirs=affichageCases(Matrice)[4]
        Ynoirs=affichageCases(Matrice)[5]
        
        plt.plot(Xverts,Yverts,"o",color="green")
        plt.plot(Xrouges,Yrouges,"o",color="red")
        plt.plot(Xnoirs,Ynoirs,"o",color="black")
        #pause de 50 ms entre deux itérations (uniquement pour l'affichage, durée pas représentative de la réalité)
        plt.pause(0.05)

plt.show()
