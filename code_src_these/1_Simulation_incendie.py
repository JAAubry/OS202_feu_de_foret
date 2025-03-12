# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 14:47:34 2021

@author: elise.foulatier
"""

import random
import numpy as np
import matplotlib.pyplot as plt

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


#############################################################################
### SIMULATION DE L'INCENDIE
#############################################################################

### DEFINITION DE L'ESPACE
#Longueur totale d'un côté de l'espace carré modélisé (en m)
L = 1000
#Espace entre deux points d'observation (en m)
longueur = 50 
#nombre de points observés sur une ligne du quadrillage
N = int(L/longueur) + 1 
#représentation du quadrillage par une matrice carrée de taille N*N (modèle bidimensionnel)
Matrice = np.zeros((N,N))

### PROBABILITES PERMETTANT LA SIMULATION DE L'INCENDIE
#probabilité que le feu se propage à un point voisin (VERT vers ROUGE)
p1=0.5
#probabilité qu'une case passe de l'état brûlant à l'état brûlé (ROUGE vers NOIR)
p2=0.3


### AFFICHAGE DU MAILLAGE AVANT DEPART DU FEU (que des points sains)
#Construction de la liste des abscisses
X=[]
for i in range(N):
    X.append(longueur*i)
    
#Affichage du maillage ligne par ligne avant départ du feu
for i in range(N):
    plt.plot(X,[i*longueur]*N,"o",color="green")
    

#Position initiale du feu (les coordonnées du point doivent être multiples de longueur = 50)
#coordonnées en m (position réelle du point de départ du feu)
x0 = 750
y0= 450
#Couple d'entiers désignant la position du point dans Matrice
m0=int(x0/longueur)
n0=int(y0/longueur)
#le point de départ est rouge donc on affecte la valeur 1 au coefficient de la matrice correspondant
Matrice[n0][m0]=1
#affichage du point de départ du feu
plt.plot([x0],[y0],"o",color="red")
#pause d'une seconde pour avoir le temps d'identifier sur la figure la position du point de départ
plt.pause(1)


### BOUCLE EFFECTUANT LA SIMULATION  DE LA PROPAGATION DU FEU TANT QUE DES CASES ROUGES SUBSISTENT


while casesRouges(Matrice)!=[]:
    #à chaque itération on répertorie les coordonnées des points rouges dans une liste
    ListeCasesRouges = casesRouges(Matrice)
    #on parcourt la liste des points rouges et on envisage pour chacun d'eux la propagation du feu dans 4 directions
    for caseRouge in ListeCasesRouges :
        #coordonnées du point brulant considéré
        #entier désignant l'abscisse
        m=caseRouge[0]
        #entier désignant l'ordonnée
        n=caseRouge[1]
        
        if m<N-1 :
            #propagation du feu dans la direction x croissants
            nombre_aleat_case_voisine=random.randint(1,100)
            #si l'entier généré est inférieur à 100 fois p1 alors la case devient rouge
            if nombre_aleat_case_voisine<p1*100 :
                if Matrice[n][m+1]==0:
                    Matrice[n][m+1]=1
        
        if m>0:
            #propagation du feu dans la direction x décroissants
            nombre_aleat_case_voisine=random.randint(1,100)
            if nombre_aleat_case_voisine<p1*100 :
                if Matrice[n][m-1]==0:
                    Matrice[n][m-1]=1
                
        if n<N-1 :
            #propagation du feu dans la direction y croissants
            nombre_aleat_case_voisine=random.randint(1,100)
            if nombre_aleat_case_voisine<p1*100 :
                if Matrice[n+1][m]==0:
                    Matrice[n+1][m]=1
        
        if n>0:
            #propagation du feu dans la direction y décroissants
            nombre_aleat_case_voisine=random.randint(1,100)
            if nombre_aleat_case_voisine<p1*100 : 
                if Matrice[n-1][m]==0:
                    Matrice[n-1][m]=1       
                        
        #évolution de l'état d'une case rouge
        nombre_aleat_reste_brulant = random.randint(1,100)
        if nombre_aleat_reste_brulant < p2*100: #si l'entier généré est inférieur à 100 fois p2 alors la case devient noire
              Matrice[n][m]=2                     
        
    #AFFICHAGE DES POINTS DANS MATPLOTLIB A CHAQUE ITERATION POUR SUIVRE L'EVOLUTION DU FEU
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


