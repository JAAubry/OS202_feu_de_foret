# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 09:39:20 2021

@author: elise.foulatier
"""

### DANS CE FICHIER ON TRACE LES COURBES DONNANT PLUSIEURS ESTIMEES DE P1 POUR PLUSIEURS VITESSES DE VENT
### LE VENT EST TOUJOURS ORIENTE D'OUEST EN EST

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
### SIMULATION DE PLUSIEURS INCENDIES
#############################################################################

### GEOMETRIE DU MODELE
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

   
### POSITION INITIALE DU FEU (coordonnées multiples de longueur = 50)
x0 = 0
y0= 500
m0=int(x0/longueur)
n0=int(y0/longueur)
Matrice[n0][m0]=1

#liste des estimations de p1 au sens du maximum de vraisemblance
listeEstimationP1=[]

### PARAMETRES LIES A LA VITESSE ET DIRECTION DU VENT
#Introduction du paramètre vmax : vitesse maximale du vent qu'on puisse rencontrer 
#elle permet de déterminer la valeur des facteurs de propagation par loi proportionnelle
Vmax = 60               #vitesse maximale en km/h

#construction listeVxVent qui contient toutes les vitesses pour lesquelles on relève plusieurs valeurs de p1
listeVxVent=[]
Vx=0                        #vitesse initiale
nombreVitesses=int(Vmax/5)  #pas de 5km/h entre 2 vitesses
for i in range (nombreVitesses+1):
    listeVxVent.append(Vx+i*5)
    
    
    
### ON PROCEDE A PLUSIEURS SIMULATIONS POUR DIFFERENTES VITESSES DU VENT
### IL Y A 5 SIMULATIONS POUR CHAQUE VITESSE
for Vx in listeVxVent :
    
    #calcul des 4 facteurs de propagation pour chaque vx
    facteur_propagation_vx_plus = (Vx/Vmax)*((1-p1)/p1) + 1
    facteur_propagation_vx_moins = (1 - Vx/Vmax)*(1/p1)

    facteur_propagation_vy_plus = 1
    facteur_propagation_vy_moins = 1
    
    l = 0
    #on effectue cinq fois la boucle pour avoir 5 estimées pour chaque Vx
    while l!=5:
        Matrice = np.zeros((N,N))
        x0 = 0
        y0= 500
        m0=int(x0/longueur)
        n0=int(y0/longueur)
        Matrice[n0][m0]=1
        
        listeProba1num=[]    #liste qui stocke les estimations du numérateur de p1 à chaque itération
        listeProba1den=[]    #liste qui stocke les estimations du numérateur de p1 à chaque itération
        k=0 #compteur des iterations
        j=False #compteur pour savoir si on a atteint l'extremite de l'espace
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
                
                for caseRouge in ListeCasesRouges :
                    m=caseRouge[0]
                    n=caseRouge[1]
                    if n==N-1:    #si le bord est du quadrillage est atteint j=Vrai et la simulation de la propagation s'arrête
                        j=True
            
        
                    if m<N-1 and j==False :
                        if Matrice[n][m+1]==0:
                            #incrémentation du nombre de points verts atteignables
                            nbCasesVertsAtteignables += 1
                            nombre_aleat_case_voisine=random.randint(1,100)
                            if nombre_aleat_case_voisine<facteur_propagation_vx_plus*p1*100 :
                                Matrice[n][m+1]=1
        
                    if m>0 and j==False:
                        if Matrice[n][m-1]==0:
                            nbCasesVertsAtteignables += 1
                            nombre_aleat_case_voisine=random.randint(1,100)
                            if nombre_aleat_case_voisine<facteur_propagation_vx_moins*p1*100 :
                                Matrice[n][m-1]=1
                
                    if n<N-1 and j==False :
                        if Matrice[n+1][m]==0:
                            nbCasesVertsAtteignables += 1
                            nombre_aleat_case_voisine=random.randint(1,100)
                            if nombre_aleat_case_voisine<facteur_propagation_vy_plus*p1*100 :
                                Matrice[n+1][m]=1
        
                    if n>0 and j==False:
                        if Matrice[n-1][m]==0:
                            nbCasesVertsAtteignables += 1
                            nombre_aleat_case_voisine=random.randint(1,100)
                            if nombre_aleat_case_voisine<facteur_propagation_vy_moins*p1*100 : 
                                Matrice[n-1][m]=1      
                        
                    #évolution de l'état d'une case rouge
                    nombre_aleat_reste_brulant = random.randint(1,100)
                    if nombre_aleat_reste_brulant < p2*100: 
                        Matrice[n][m]=2                     

                anciensCasesRouges=nbCasesRougesInit - (compterCasesNoires(Matrice)-nbCasesNoirsInit)
                nouveauxCasesRouges=len(casesRouges(Matrice)) - anciensCasesRouges
                k+=1
                
                #ajout de nouvelles valeurs aux sommes au numérateur et au dénominateur de p1
                if nbCasesVertsInit>0:
                    listeProba1num.append(nouveauxCasesRouges)
                    listeProba1den.append(nbCasesVertsAtteignables)
        
        #CALCUL ESTIMATION P1 A CHAQUE ITERATION
        estimation_p1 = sum(listeProba1num)/sum(listeProba1den)
        listeEstimationP1.append(estimation_p1)
        #une des cinq simulations a été effectuée donc on incrémente l de 1
        l+=1
        
        
        
###CONSTRUCTION DE LA LISTE DE VITESSES POUR L'AFFICHAGE DANS MATPLOTLIB     
listeVitesseX=[]
for Vx in listeVxVent :
    #on ajoute 5 fois Vx comme il y a 5 mesures de p1 par vitesse
    listeVitesseX.append(Vx)
    listeVitesseX.append(Vx)
    listeVitesseX.append(Vx)
    listeVitesseX.append(Vx)
    listeVitesseX.append(Vx)
        
    
### AFFICHAGE DE LA COURBE DONNANT LES ESTIMEES DE P1 EN FONCTION DE LA VITESSE DU VENT
plt.plot(listeVitesseX,listeEstimationP1,"+",color="red")
plt.xlabel("Vitesse du vent Vx (m/s)")
plt.ylabel("Estimées de p1")
plt.grid()
plt.show()