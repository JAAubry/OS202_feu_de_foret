# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 09:38:16 2021

@author: elise.foulatier
"""

### DANS CE FICHIER ON MODIFIE LA MANIERE DONT SE PROPAGE LE FEU AVEC UN VENT ALLANT D'OUEST EN EST
### ON REALISE UNE ESTIMATION AU SENS DU MAXIMUM DE VRAISEMBLANCE DE P1 EN FIN DE SIMULATION

import random
import numpy as np
import matplotlib.pyplot as plt

###################################################################################################
### DEFINITION D'UN CERTAIN NOMBRE DE FONCTIONS UTILES POUR LA SIMULATION ET L'ESTIMATION DE P1
###################################################################################################

#Fonction prenant en paramètre une matrice et retournant toutes les coordonnées des points rouges
def casesRouges(Matrice):
    M=len(Matrice[0]) #nombre abscisses
    N=len(Matrice) #nombre d'ordonnees
    
    ListePointsBrulants=[] #liste qui contient les points brulants
    
    for i in range(N):
        for j in range(M):
            if Matrice[i][j]==1:
                ListePointsBrulants.append([j,i])
    return ListePointsBrulants

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


#Fonction qui prend en entrée la matrice modélisant l'espace et retourne le nombre de cases vertes qu'elle contient
def compterCasesVertes(Matrice):
    M=len(Matrice[0]) #nombre abscisses
    N=len(Matrice) #nombre d'ordonnees
    
    sommeCasesVertes=0
    #balayage de la matrice
    for i in range(N):
        for j in range(M):
            if Matrice[i][j]==0: #cases vertes représentées par un 0
                sommeCasesVertes+=1
    return sommeCasesVertes


#Fonction qui prend en entrée la matrice modélisant l'espace et retourne le nombre de cases vertes qu'elle contient
def compterCasesNoires(Matrice):
    M=len(Matrice[0]) #nombre abscisses
    N=len(Matrice) #nombre d'ordonnees
    
    sommePointsBrules=0 #
    
    for i in range(N):
        for j in range(M):
            if Matrice[i][j]==2:     #cases vertes représentées par un 2
                sommePointsBrules+=1
    return sommePointsBrules


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
p1=0.4
#probabilité qu'une case passe de l'état brûlant à l'état brûlé (ROUGE vers NOIR)
p2=0.3


### PARAMETRES LIES A LA VITESSE ET DIRECTION DU VENT
#Introduction du paramètre vmax : vitesse maximale du vent qu'on puisse rencontrer 
#elle permet de déterminer la valeur des facteurs de propagation par loi proportionnelle
Vmax = 60        #correspond à un vent de 60 km/h

#ici vent de direction ouest vers est
Vx = 45     #valeur choisie arbitrairement pour cette simulation : peut être modifiée 
Vy = 0

#Calcul des 4 facteurs d'après la loi décrit dans le rapport
facteur_propagation_vx_plus = (Vx/Vmax)*((1-p1)/p1) + 1
facteur_propagation_vx_moins = (1 - Vx/Vmax)*(1/p1)
#Le feu ne se propage pas différemment latéralement donc on ne modifie pas la valeur de la probabilité dans la direction y
facteur_propagation_vy_plus = 1
facteur_propagation_vy_moins = 1


### AFFICHAGE DU MAILLAGE AVANT DEPART DU FEU (que des points verts)
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

# listes qui stockent à chaque itération le k-ième terme de la somme au numérateur et de la somme au dénominateur
# de l'estimée au maximum de vraisemblance de p1
listeProba1num=[]
listeProba1den=[]
k=0 #compteur des itérations de la boucle
j=False #variable booleenne qui permet de savoir si on a atteint l'extrêmité de l'espace


while casesRouges(Matrice)!=[]:
    ListeCasesRouges = casesRouges(Matrice)
    #variable donnant le nombre de cases rouges à l'itération k
    nbCasesRougesInit=len(ListeCasesRouges)
    #variable donnant le nombre de cases noirs à l'itération k
    nbCasesNoirsInit=compterCasesNoires(Matrice)
    #variable donnant le nombre de cases verts à l'itération k : 
    #dernier état possible donc on soustrait au nombre de cases total le nombre de cases rouges ou noires
    nbCasesVertsInit = N**2 - (nbCasesRougesInit+nbCasesNoirsInit)
    #intialisation d'une variable donnant le nombre de points verts qu'on peut atteindre en une itération 
    #(points voisins d'un point rouge et qui n'est pas déjà rouge)
    nbCasesVertsAtteignables = 0
    #affichage du nombre de points rouges
    print("nombre Cases rouges :", nbCasesRougesInit)
    #on parcourt la liste des points brulants et on envisage pour chacun d'eux la propagation du feu dans 4 directions
    for caseRouge in ListeCasesRouges :
        #coordonnées du point rouge considéré
        m=caseRouge[0]
        n=caseRouge[1]
        if m==N-1: #si le bord est du quadrillage est atteint j=Vrai et la simulation de la propagation s'arrête
            j=True
            
        
        if m<N-1 and j==False :
            if Matrice[n][m+1]==0:
                #incrémentation du nombre de points verts atteignables
                nbCasesVertsAtteignables += 1
                #propagation du feu dans la direction x croissants
                nombre_aleat_case_voisine=random.randint(1,100)
                if nombre_aleat_case_voisine<facteur_propagation_vx_plus*p1*100 :
                    Matrice[n][m+1]=1
        
        if m>0 and j==False:
            if Matrice[n][m-1]==0:
                nbCasesVertsAtteignables += 1
                #propagation du feu dans la direction x décroissants
                nombre_aleat_case_voisine=random.randint(1,100)
                if nombre_aleat_case_voisine<facteur_propagation_vx_moins*p1*100 :
                    Matrice[n][m-1]=1
                
        if n<N-1 and j==False :
            if Matrice[n+1][m]==0:
                nbCasesVertsAtteignables += 1
                #propagation du feu dans la direction y croissants
                nombre_aleat_case_voisine=random.randint(1,100)
                if nombre_aleat_case_voisine<facteur_propagation_vy_plus*p1*100 :
                    Matrice[n+1][m]=1
        
        if n>0 and j==False:
            if Matrice[n-1][m]==0:
                nbCasesVertsAtteignables += 1
                #propagation du feu dans la direction y décroissants
                nombre_aleat_case_voisine=random.randint(1,100)
                if nombre_aleat_case_voisine<facteur_propagation_vy_moins*p1*100 : 
                    Matrice[n-1][m]=1      
                        
        #evolution de l'etat d'une case rouge
        nombre_aleat_reste_brulant = random.randint(1,100)
        if nombre_aleat_reste_brulant < p2*100: #si le nombre généré est inférieur à 100 fois la probabilité alors la case passe à l'état brûlé
              Matrice[n][m]=2                     #la case considérée est désormais brûlée, coef égal à 2

    #variable donnant le nombre de cases passées de rouge à noir
    anciensCasesRouges=nbCasesRougesInit - (compterCasesNoires(Matrice)-nbCasesNoirsInit)
    #variable donnant le nombre de cases passées de vert à rouge
    nouveauxCasesRouges=len(casesRouges(Matrice)) - anciensCasesRouges
    
    #on augmente l'instant de 1 (on est à temps discret)
    k+=1
    #affichage dans la console des observations
    print("Nombre de cases vertes à l'itération k=",k,":",compterCasesVertes(Matrice))
    print("Nombre de cases vertes qui étaient atteignables par le feu avant l'itération k=",k,":",nbCasesVertsAtteignables)
    print("Nombre de cases restées rouges à l'itération k =",k,":",anciensCasesRouges)
    print("Nombre de nouvelles cases rouges à l'itération k=",k,":",nouveauxCasesRouges)
    print("Nombre de cases noires à l'itération k =",k,":",compterCasesNoires(Matrice))
        
    #ESTIMATION P1 AVEC VENT
    if nbCasesVertsInit>0:
        listeProba1num.append(nouveauxCasesRouges)
        listeProba1den.append(nbCasesVertsAtteignables)
        
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

#calcul de l'estimation de p1 au sens du maximum de vraisemblance
#à l'aide de l'expression donnée dans le rapport
estimation_p1 = sum(listeProba1num)/sum(listeProba1den)
print("Estimation de p1 au maximum de vraisemblance :   ",estimation_p1)