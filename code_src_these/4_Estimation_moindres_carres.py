# -*- coding: utf-8 -*-
"""
Created on Tue May 11 17:33:07 2021

@author: elise.foulatier
"""

### DANS CE FICHIER ON RECUPERE LA LISTE DES MESURES DE P1
###ON CHERCHE A TROUVER UN MODELE PARAMETRIQUE SIMPLE DE P1 = F(V)

from math import log, exp
import numpy as np
import matplotlib.pyplot as plt

listeEstimationP1 = [0.16666666666666666,
 0.4382940108892922,
 0.442794498704405,
 0.44068043035766213,
 0.4497628763648395,
 0.45052304212609556,
 0.4602721685689201,
 0.45772796457727966,
 0.4684390800848149,
 0.4601822672685461,
 0.46855476884631403,
 0.45679558011049726,
 0.4636622005684125,
 0.4677999751830252,
 0.4678494623655914,
 0.48177339901477834,
 0.48534454153301626,
 0.47222498260955975,
 0.47320422174403115,
 0.4934684920054342,
 0.5013629283489096,
 0.0,
 0.48766564729867484,
 0.49347316746625014,
 0.4979186679474864,
 0.5069410258934725,
 0.49842359050445106,
 0.5034336009285231,
 0.5068570143884892,
 0.5115369181380417,
 0.5210278028706213,
 0.5159710174507602,
 0.5164948453608248,
 0.5089121745705952,
 0.5160876431109357,
 0.5256446585434967,
 0.5336194563662375,
 0.5411553110452815,
 0.5299678915447734,
 0.5234283483014448,
 0.5426160337552742,
 0.5610453246222948,
 0.5402127262548984,
 0.5537465486051605,
 0.5394137452960982,
 0.5581765173601908,
 0.16666666666666666,
 0.5574386133846895,
 0.5782530056388978,
 0.5649429741999779,
 0.5687344438062416,
 0.5686968838526912,
 0.5709683717297931,
 0.5721233253808038,
 0.5627927146574154,
 0.6005920078934386,
 0.6307270591306617,
 0.6017146486028789,
 0.5927367055771725,
 0.6235309652713588,
 0.6695869837296621,
 0.6426574357608841,
 0.6512820512820513,
 0.6751427088738973,
 0.6624671189256541]



##########################################################################
### MODELE EXPONENTIEL
### P1 = lambda.exp(alpha.v)
##########################################################################

#construction du vecteur y tilde : on applique ln à chaque mesure (chaque élément de listeEstimationP1)
listeEstimationP1ln=[]
for estimation in listeEstimationP1 :
    if estimation!=0:
        listeEstimationP1ln.append(log(estimation))
    #si des mesures sont égales à 0 on ne peut pas appliquer ln donc on leur affecte la valeur 0
    else :
        listeEstimationP1ln.append(0)
        
        
### SANS PONDERATION

#Définition de la matrice R des régresseurs
R=[] 
nbMesuresParVitesse=5      #on effectue 5 mesures pour chaque vitesse donnée
intervalleVitesse=5        #intervalle entre deux vitesses : ici 5 km/h
j=0
while j<=60:
    l=0
    while l<nbMesuresParVitesse:
        R.append([1,j]) 
        l+=1
    j+=intervalleVitesse
        
#Calcul matriciel pour estimer q au sens des moindres carres : voir rapport

#matrice pseudo-inverse
matPseudoInv=np.linalg.inv(np.dot(np.transpose(R),R))
#second membre en tR.y
secondMembre = np.dot(np.transpose(R),listeEstimationP1ln)
#estimee qMC du vecteur des paramètres
qMC=np.dot(matPseudoInv,secondMembre)

#approximation de l'estimee de p au sens des moindres carrés
pMC=[exp(qMC[0]),qMC[1]]


### AFFICHAGE DES MESURES ET DU MODELE ESTIME DANS MATPLOTLIB
#construction de la liste des abscisses pour l'affichage des mesures
listeVxVentMesures=[]

Vmax=60
nbVitesses=int(Vmax/intervalleVitesse)
for i in range (nbVitesses+1):
    listeVxVentMesures.append(i*intervalleVitesse)
    listeVxVentMesures.append(i*intervalleVitesse)
    listeVxVentMesures.append(i*intervalleVitesse)
    listeVxVentMesures.append(i*intervalleVitesse)
    listeVxVentMesures.append(i*intervalleVitesse)
    
#construction de la liste des abscisses pour l'affichage du modèle
listeVxVentModele=[]
nbAbscissesModele=Vmax*10
for i in range(nbAbscissesModele):
    listeVxVentModele.append(i*0.1)
#construction des valeurs prises par le modèle
p1Modele=[]
for vitesse in listeVxVentModele:
    p1Modele.append(pMC[0]*exp(vitesse*pMC[1]))
    
    
#Affichage des courbes  
plt.plot(listeVxVentMesures,listeEstimationP1,"+",color="red")
plt.plot(listeVxVentModele,p1Modele,color="green")
plt.xlabel("Vitesse du vent Vx (m/s)")
plt.ylabel("Estimées de p1")
plt.grid()
plt.show()

#Calcul du critère : somme des carrés de l'erreur
n=len(listeEstimationP1)
critere=0
for i in range(n):
    critere+=(listeEstimationP1[i]-pMC[0]*exp(listeVxVentMesures[i]*pMC[1]))**2
print("Valeur du critère =",critere)


### AVEC PONDERATION

#matrice R
R=[]
nbMesuresParVitesse=5
n=len(listeEstimationP1)
j=0
while j<=60:
    l=0
    while l<nbMesuresParVitesse:
        R.append([1,j]) 
        l+=1
    j+=nbMesuresParVitesse
    
#Définition de la matrice W de pondération
#W est une matrice diagonale contenant des 1 sur sa diagonale
#si une des mesures est aberrante (trop éloignée des autres) : terme mis à 0.001

diagonaleW=[1]*(nbMesuresParVitesse*(nbVitesses+1))
W=np.diag(diagonaleW)
n=len(listeEstimationP1)

for i in range(n-1):
    if abs(listeEstimationP1[i])<0.4: #mesures aberrantes quand inférieures à 0.4 d'après le graphique des mesures
                                      #et comme on a fixé p1 = 0.4, on ne peut pas avoir des mesures plus faibles
        W[i][i]=0.001
        
        
#Calcul matriciel pour estimer q au sens des moindres carrés : voir rapport

#matrice à inverser
matriceInversable = np.dot(np.dot(np.transpose(R),W),R)
matriceInversee = np.linalg.inv(matriceInversable)
#definition du second facteur dans l'expression de pMC
secondMembre = np.dot(np.dot(np.transpose(R),W),listeEstimationP1ln)
#estimee de q au sens des MC
qMC = np.dot(matriceInversee,secondMembre)

#approximation de l'estimée de p au sens des moindres carrés
pMC=[exp(qMC[0]),qMC[1]]


### AFFICHAGE DES MESURES ET DU MODELE ESTIME DANS MATPLOTLIB
#construction de la liste des abscisses pour l'affichage des mesures
listeVxVentMesures=[]

Vmax=60
nbVitesses=int(Vmax/intervalleVitesse)
for i in range (nbVitesses+1):
    listeVxVentMesures.append(i*intervalleVitesse)
    listeVxVentMesures.append(i*intervalleVitesse)
    listeVxVentMesures.append(i*intervalleVitesse)
    listeVxVentMesures.append(i*intervalleVitesse)
    listeVxVentMesures.append(i*intervalleVitesse)
    
#construction de la liste des abscisses pour l'affichage du modèle
listeVxVentModele=[]
nbAbscissesModele=Vmax*10
for i in range(nbAbscissesModele):
    listeVxVentModele.append(i*0.1)
#construction des valeurs prises par le modèle
p1Modele=[]
for vitesse in listeVxVentModele:
    p1Modele.append(pMC[0]*exp(vitesse*pMC[1]))
    
    
#Affichage des courbes  
plt.plot(listeVxVentMesures,listeEstimationP1,"+",color="red")
plt.plot(listeVxVentModele,p1Modele,color="green")
plt.xlabel("Vitesse du vent Vx (m/s)")
plt.ylabel("Estimées de p1")
plt.grid()
plt.show()


#Calcul du critère : somme des carrés de l'erreur
critere=0
for i in range(n):
    critere+=W[i][i]*(listeEstimationP1[i]-pMC[0]*exp(listeVxVentMesures[i]*pMC[1]))**2
print("Valeur du critère =",critere)




##########################################################################
### MODELE POLYNOMIAL
### P1 = p0 + p1v+...+p4v^4
##########################################################################



### SANS PONDERATION

#définition de la matrice R des régresseurs pour un polynôme d'ordre 3
R=[]
nbMesuresParVitesse=5
n=len(listeEstimationP1)
j=0
while j<=nbVitesses:
    v=listeVxVentMesures[j*5]
    l=0
    while l<nbMesuresParVitesse:
        R.append([1,v,v**2,v**3]) 
        l+=1
    j+=1

        
#Calcul matriciel pour estimer p au sens des moindres carrés : voir rapport

#matrice pseudo-inverse
matPseudoInv=np.linalg.inv(np.dot(np.transpose(R),R))
#second membre en tR.y
secondMembre = np.dot(np.transpose(R),listeEstimationP1)
#estimee pMC du vecteur des paramètres
pMC=np.dot(matPseudoInv,secondMembre)


### AFFICHAGE DES MESURES ET DU MODELE ESTIME DANS MATPLOTLIB
listeVxVentModele=[]
nbAbscissesModele=Vmax*10
for i in range(nbAbscissesModele):
    listeVxVentModele.append(i*0.1)
#construction du modèle
p1Modele=[]
for vitesse in listeVxVentModele:
    p1Modele.append(pMC[0] + pMC[1]*vitesse +pMC[2]*(vitesse**2)+pMC[3]*(vitesse**3))
    
    
#Affichage des courbes  
plt.plot(listeVxVentMesures,listeEstimationP1,"+",color="red")
plt.plot(listeVxVentModele,p1Modele,color="green")
plt.xlabel("Vitesse du vent Vx (m/s)")
plt.ylabel("Estimées de p1")
plt.grid()
plt.show()

#Calcul du critère des moindres carrés
n=len(listeEstimationP1)
critere=0
for i in range(n):
    critere+=(listeEstimationP1[i]-(pMC[0] + pMC[1]*listeVxVentMesures[i]+pMC[2]*(listeVxVentMesures[i]**2)+pMC[3]*(listeVxVentMesures[i]**3)))**2
print("Valeur du critère = ",critere)



### AVEC PONDERATION


#définition de la matrice R des régresseurs pour un polynôme d'ordre 3
R=[]
nbMesuresParVitesse=5
n=len(listeEstimationP1)
j=0
while j<=nbVitesses:
    v=listeVxVentMesures[j*5]
    l=0
    while l<nbMesuresParVitesse:
        R.append([1,v,v**2,v**3]) 
        l+=1
    j+=1
    
#definition de la matrice W de ponderation
#W est une matrice diagonale contenant des 1 sur sa diagonale
#si une des mesures est aberrante (trop éloignée des autres) : terme mis à 0.001

diagonaleW=[1]*(nbMesuresParVitesse*(nbVitesses+1))
W=np.diag(diagonaleW)
n=len(listeEstimationP1)
for i in range(n-1):
    if abs(listeEstimationP1[i])<0.4:
        W[i][i]=0.001


#Calcul matriciel pour estimer p au sens des moindres carrés : voir rapport

#matrice à inverser
matriceInversable = np.dot(np.dot(np.transpose(R),W),R)
matriceInversee = np.linalg.inv(matriceInversable)
#definition du second facteur dans l'expression de pMC
secondMembre = np.dot(np.dot(np.transpose(R),W),listeEstimationP1)
#estimee de p au sens des MC
pMC = np.dot(matriceInversee,secondMembre)


### AFFICHAGE DES MESURES ET DU MODELE ESTIME DANS MATPLOTLIB

listeVxVentModele=[]
nbAbscissesModele=Vmax*10
for i in range(nbAbscissesModele):
    listeVxVentModele.append(i*0.1)
#construction du modèle
p1Modele=[]
for vitesse in listeVxVentModele:
    p1Modele.append(pMC[0] + pMC[1]*vitesse+pMC[2]*(vitesse**2)+pMC[3]*(vitesse**3))
    
    
#Affichage des courbes  
plt.plot(listeVxVentMesures,listeEstimationP1,"+",color="red")
plt.plot(listeVxVentModele,p1Modele,color="green")
plt.xlabel("Vitesse du vent Vx (m/s)")
plt.ylabel("Estimées de p1")
plt.grid()
plt.show()

#calcul de la somme des carrés de l'erreur
n=len(listeEstimationP1)
critere=0
for i in range(n):
    critere+=W[i][i]*(listeEstimationP1[i]-(pMC[0] + pMC[1]*listeVxVentMesures[i]+pMC[2]*(listeVxVentMesures[i]**2)+pMC[3]*(listeVxVentMesures[i]**3)))**2
print("Valeur du critère =",critere)