#install.packages('tidyverse')
#install.packages("cowplot")
#install.packages("gridExtra")
library('tidyverse')
library("cowplot")
library("gridExtra")


#### Fonctions ####
MAR = function(P, alpha, beta, N, debut){
  Mar = debut
  for(i in 1:(N-1)){
    if((P[Mar[i],1] > 10**300) | (P[Mar[i],2] > 10**300)){
      warning('renforcement trop important')
      return(list(Mar = Mar, P = P))
      }
    
    #On tire un nombre aleatoire entre 0 et 1
    U = runif(1)
    
    #Si on est au bord de gauche et que notre U vaut 1
    #La boucle d'apres voudra aller a gauche (impossible)
    # Il faut donc s'assurer ici que nous ne sommes pas dans ce cas
    if((U == 1) && (Mar[i] == 1)){                                
      P[1,2] = beta*P[1,2] + alpha                           
      Mar = c(Mar, Mar[i] + 1)
    }
    else if(U < (P[Mar[i],2] / sum(P[Mar[i],]))){                 
      P[Mar[i],2] = beta*P[Mar[i],2] + alpha                        #On va a droite et on ajoute delta 
      Mar = c(Mar, Mar[i] + 1)
    }
    else{                                                                       #Sinon pareil mais on va a gauche
      P[Mar[i],1] = beta*P[Mar[i],1] + alpha
      Mar = c(Mar, Mar[i] - 1)
    }
  }
  return(list(Mar = Mar, P = P))
}


MAR2d = function(P, Q, alpha, beta, N, debut){
  Mar = matrix(0, nrow = 2, ncol = N)
  Mar[1,1] = debut
  Mar[2,1] = debut
  K = nrow(P)
  Sommet = matrix(0, nrow = K, ncol = K)                         # On compte chaque passage à chaque sommet
  ArreteP = matrix(0, nrow = K, ncol = 2*K)                      # On compte chaque passage à chaque arête horizontale
  ArreteQ = matrix(0, nrow = K, ncol = 2*K)                      # On compte chaque passage à chaque arête verticale
  for(i in 1:(N-1)){
    U = runif(1)                                                  #On tire un nombre aleatoire entre 0 et 1
    a = P[Mar[1,i],  Mar[2,i]]                                       #On créer des variables pour simplifier l'écriture
    b = P[Mar[1,i],  Mar[2,i] + K]
    c = Q[Mar[1,i],  Mar[2,i] + K]
    d = Q[Mar[1,i],  Mar[2,i]]
    Sommet[Mar[1,i],  Mar[2,i]] = Sommet[Mar[1,i],  Mar[2,i]] + 1
    if((a > 10**300) | (b > 10**300) | (c > 10**300) | (d > 10**300)){
      warning('renforcement trop important')
      return(list(Mar = matrix(Mar[Mar > 0], nrow = 2, byrow = FALSE), P = P, Q = Q, i = i, ArreteP = ArreteP, ArreteQ = ArreteQ, Sommet = Sommet))
    }
    if(U < (a / (a + b + c + d))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
      P[Mar[1,i], Mar[2,i]] = beta * P[Mar[1,i], Mar[2,i]] + alpha                           #On va a gauche et on ajoute delta 
      Mar[1, i + 1] = Mar[1, i] - 1
      Mar[2, i + 1] = Mar[2, i]
      ArreteP[Mar[1,i], Mar[2,i]] = ArreteP[Mar[1,i], Mar[2,i]] + 1
    }
    else if(U < ((a + b) / (a + b + c + d))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
      P[Mar[1,i], Mar[2,i] + K] = beta * P[Mar[1,i], Mar[2,i] + K] + alpha                                               #On va a droite et on enleve delta
      Mar[1, i + 1] = Mar[1, i] + 1
      Mar[2, i + 1] = Mar[2, i]
      ArreteP[Mar[1,i], Mar[2,i] + K] = ArreteP[Mar[1,i], Mar[2,i] + K] + 1
    }
    else if(U < ((a + b + c) / (a + b + c + d))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
      Q[Mar[1,i], Mar[2,i] + K] = beta * Q[Mar[1,i], Mar[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
      Mar[1, i + 1] = Mar[1, i]
      Mar[2, i + 1] = Mar[2, i] + 1
      ArreteQ[Mar[1,i], Mar[2,i] + K] = ArreteQ[Mar[1,i], Mar[2,i] + K] + 1
    }
    else if((U != 1) | (Mar[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
      Q[Mar[1,i], Mar[2,i]] = beta * Q[Mar[1,i], Mar[2,i]] + alpha                                                      #On va derriere et on enleve delta
      Mar[1, i + 1] = Mar[1, i]
      Mar[2, i + 1] = Mar[2, i] - 1
      ArreteQ[Mar[1,i], Mar[2,i]] = ArreteQ[Mar[1,i], Mar[2,i]] + 1
    }
    else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
      Mar[1, i + 1] = Mar[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
      Mar[2, i + 1] = Mar[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
    } 
  }
  Sommet[Mar[1,N],  Mar[2,N]] = Sommet[Mar[1,N],  Mar[2,N]] + 1
  return(list(Mar = Mar, P = P, Q = Q, ArreteP = ArreteP, ArreteQ = ArreteQ, Sommet = Sommet))
}


MARri = function(P, Q, alpha, beta, N, debut){
  Mar = matrix(0, nrow = 2, ncol = N)
  Mar[1,1] = debut
  Mar[2,1] = debut
  K = nrow(P)
  Sommet = matrix(0, nrow = K, ncol = K)                         # On compte chaque passage à chaque sommet
  ArreteP = matrix(0, nrow = K, ncol = 2*K)                      # On compte chaque passage à chaque arête horizontale
  ArreteQ = matrix(0, nrow = K, ncol = 2*K)                      # On compte chaque passage à chaque arête verticale
  a = P[debut,debut]                                                  #On initialise 4 variable a, b, c et d qui contiennent
  b = P[debut,debut + K]                                                  # le poids d'aller dans chacune des directions par rapport
  c = Q[debut,debut + K]                                                  # a l'endroit où nous sommes a l'instant t
  d = Q[debut,debut]                                                  
  e = 0                                                           #On initialise un e qui contiendra le mouvement realise precedemment (1 si on est alle a gauche, 2 a droite, 3 devant et 4 derriere)
  for(i in 1:(N-1)){
    if((P[Mar[1,i], Mar[2,i]] > 10**300) | (P[Mar[1,i], Mar[2,i] + K] > 10**300) | (Q[Mar[1,i], Mar[2,i]] > 10**300) | (Q[Mar[1,i], Mar[2,i] + K] > 10**300)){
      warning('renforcement trop important')
      return(list(Mar = matrix(Mar[Mar > 0], nrow = 2, byrow = FALSE), P = P, Q = Q, i = i, ArreteP = ArreteP, ArreteQ = ArreteQ, Sommet = Sommet))
    }
    Sommet[Mar[1,i],  Mar[2,i]] = Sommet[Mar[1,i],  Mar[2,i]] + 1
    U = runif(1)                                                  #On tire un nombre aleatoire entre 0 et 1
    if(U < (a / (a + b + c + d))){                                #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
      P[Mar[1,i], Mar[2,i]] = beta * P[Mar[1,i], Mar[2,i]] + alpha                       #On va a gauche et on ajoute delta 
      Mar[1, i + 1] = Mar[1, i] - 1
      Mar[2, i + 1] = Mar[2, i]
      ArreteP[Mar[1,i], Mar[2,i]] = ArreteP[Mar[1,i], Mar[2,i]] + 1
      e = 1                                                       #On enregistre notre mouvement (on est alle a gauche)
    }
    else if(U < ((a + b) / (a + b + c + d))){                     #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
      P[Mar[1,i], Mar[2,i] + K] = beta * P[Mar[1,i], Mar[2,i] + K] + alpha                       #On va a droite et on enleve delta
      Mar[1, i + 1] = Mar[1, i] + 1
      Mar[2, i + 1] = Mar[2, i]
      ArreteP[Mar[1,i], Mar[2,i] + K] = ArreteP[Mar[1,i], Mar[2,i] + K] + 1
      e = 2                                                       #On enregistre notre mouvement (on est alle a droite)
    }
    else if(U < ((a + b + c) / (a + b + c + d))){                 #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
      Q[Mar[1,i], Mar[2,i] + K] = beta * Q[Mar[1,i], Mar[2,i] + K] + alpha                       #On va devant et on ajoute delta
      Mar[1, i + 1] = Mar[1, i]
      Mar[2, i + 1] = Mar[2, i] + 1
      ArreteQ[Mar[1,i], Mar[2,i] + K] = ArreteQ[Mar[1,i], Mar[2,i] + K] + 1
      e = 3                                                       #On enregistre notre mouvement (on est alle devant)
    }
    else if((U != 1) | (Mar[2,i] != 1)){                          #Si U est different de 1 ou que nous ne sommes pas tout en bas
      Q[Mar[1,i], Mar[2,i]] = beta * Q[Mar[1,i], Mar[2,i]] + alpha                       #On va derriere et on enleve delta
      Mar[1, i + 1] = Mar[1, i]
      Mar[2, i + 1] = Mar[2, i] - 1
      ArreteQ[Mar[1,i], Mar[2,i]] = ArreteQ[Mar[1,i], Mar[2,i]] + 1
      e = 4                                                       #On enregistre notre mouvement (on est alle derriere)
    }
    else{                                                         #Si U est egal a 1 et que nous sommes tout en bas,                                          
      Mar[1, i + 1] = Mar[1, i]                                   #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
      Mar[2, i + 1] = Mar[2, i]                                   #Ce cas n'a quasiment aucune chance d'arrive mais il faut tout de meme l'anticiper
    } 
    a = P[Mar[1, i + 1],  Mar[2,i + 1]]                                       #On actualise nos variables avec les nouvelles valeurs pour l'instant suivant
    b = P[Mar[1, i + 1],  Mar[2,i + 1] + K]
    c = Q[Mar[1, i + 1],  Mar[2,i + 1] + K]
    d = Q[Mar[1, i + 1],  Mar[2,i + 1]]
    if(e == 1){                                                   #On met en place l'interdiction de retour en arriere
      b = 0                                                       #Si on est alle a gauche, on a l'interdiction d'aller a droite a l'instant suivant
    }
    else if(e == 2){
      a = 0                                                       #Si on est alle a droite, on a l'interdiction d'aller a gauche a l'instant suivant
    }
    else if(e == 3){
      d = 0                                                       #Si on est alle devant, on a l'interdiction d'aller derriere a l'instant suivant
    }
    else if(e ==4){
      c = 0                                                       #Si on est alle derriere, on a l'interdiction d'aller devant a l'instant suivant
    }
  }
  Sommet[Mar[1,N],  Mar[2,N]] = Sommet[Mar[1,N],  Mar[2,N]] + 1
  return(list(Mar = Mar, P = P, Q = Q, ArreteP = ArreteP, ArreteQ = ArreteQ, Sommet = Sommet))
}


MARstop = function(P, Q, alpha, beta, N, debut){
  Mar = matrix(0, nrow = 2, ncol = N)
  Mar[1,1] = debut
  Mar[2,1] = debut
  K = nrow(P)
  Sommet = matrix(0, nrow = K, ncol = K)                         # On compte chaque passage à chaque sommet
  ArreteP = matrix(0, nrow = K, ncol = 2*K)                      # On compte chaque passage à chaque arête horizontale
  ArreteQ = matrix(0, nrow = K, ncol = 2*K)                      # On compte chaque passage à chaque arête verticale
  for(i in 1:(N-1)){
    if((Mar[1,i] == K) | (Mar[1,i] == 1) | (Mar[2,i] == K) | (Mar[2,i] == 1)){                   #Si on est sur un bord
      return(list(Mar = matrix(Mar[Mar > 0], nrow = 2, byrow = FALSE), P = P, Q = Q, i = i, ArreteP = ArreteP, ArreteQ = ArreteQ, Sommet = Sommet))     #On renvoie la trajectoire
    }
    Sommet[Mar[1,i],  Mar[2,i]] = Sommet[Mar[1,i],  Mar[2,i]] + 1
    a = P[Mar[1,i],  Mar[2,i]]                                       #On créer des variables pour simplifier l'écriture
    b = P[Mar[1,i],  Mar[2,i] + K]
    c = Q[Mar[1,i],  Mar[2,i] + K]
    d = Q[Mar[1,i],  Mar[2,i]]
    if((a > 10**300) | (b > 10**300) | (c > 10**300) | (d > 10**300)){
      warning('renforcement trop important')
      return(list(Mar = matrix(Mar[Mar > 0], nrow = 2, byrow = FALSE), P = P, Q = Q, i = i, ArreteP = ArreteP, ArreteQ = ArreteQ, Sommet = Sommet))
    }
    U = runif(1)                                                      #On tire un nombre aleatoire entre 0 et 1
    if(U < (a / (a + b + c + d))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
      P[Mar[1,i], Mar[2,i]] = beta * P[Mar[1,i], Mar[2,i]] + alpha                           #On va a gauche et on ajoute delta 
      Mar[1, i + 1] = Mar[1, i] - 1
      Mar[2, i + 1] = Mar[2, i]
      ArreteP[Mar[1,i], Mar[2,i]] = ArreteP[Mar[1,i], Mar[2,i]] + 1
    }
    else if(U < ((a + b) / (a + b + c + d))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
      P[Mar[1,i], Mar[2,i] + K] = beta * P[Mar[1,i], Mar[2,i] + K] + alpha                                               #On va a droite et on enleve delta
      Mar[1, i + 1] = Mar[1, i] + 1
      Mar[2, i + 1] = Mar[2, i]
      ArreteP[Mar[1,i], Mar[2,i] + K] = ArreteP[Mar[1,i], Mar[2,i] + K] + 1
    }
    else if(U < ((a + b + c) / (a + b + c + d))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
      Q[Mar[1,i], Mar[2,i] + K] = beta * Q[Mar[1,i], Mar[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
      Mar[1, i + 1] = Mar[1, i]
      Mar[2, i + 1] = Mar[2, i] + 1
      ArreteQ[Mar[1,i], Mar[2,i] + K] = ArreteQ[Mar[1,i], Mar[2,i] + K] + 1
    }
    else if((U != 1) | (Mar[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
      Q[Mar[1,i], Mar[2,i]] = beta * Q[Mar[1,i], Mar[2,i]] + alpha                                                      #On va derriere et on enleve delta
      Mar[1, i + 1] = Mar[1, i]
      Mar[2, i + 1] = Mar[2, i] - 1
      ArreteQ[Mar[1,i], Mar[2,i]] = ArreteQ[Mar[1,i], Mar[2,i]] + 1
    }
    else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
      Mar[1, i + 1] = Mar[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
      Mar[2, i + 1] = Mar[2, i]                                 #Ce cas n'a quasiment aucune chance d'arrive mais il faut tout de meme l'anticiper
    } 
  }
  Sommet[Mar[1,N],  Mar[2,N]] = Sommet[Mar[1,N],  Mar[2,N]] + 1
  return(list(Mar = Mar, P = P,Q = Q, ArreteP = ArreteP, ArreteQ = ArreteQ, Sommet = Sommet))
}


MAR2pop = function(P, Q, R, S, alpha, beta, gamma, delta, N, xdebut, ydebut, xdebut1, ydebut1){
  Mar1 = matrix(0, nrow = 2, ncol = N)
  Mar2 = matrix(0, nrow = 2, ncol = N)
  Mar1[1,1] = xdebut
  Mar2[1,1] = xdebut1
  Mar1[2,1] = ydebut
  Mar2[2,1] = ydebut1
  K = nrow(P)
  for(i in 1:(N-1)){
    W = runif(1)
    
    U = runif(1)                                                  #On tire un nombre aleatoire entre 0 et 1
    a = P[Mar1[1,i], Mar1[2,i]]                                       #On créer des variables pour simplifier l'écriture
    b = P[Mar1[1,i], Mar1[2,i] + K]
    c = Q[Mar1[1,i], Mar1[2,i] + K]
    d = Q[Mar1[1,i], Mar1[2,i]]
    
    V = runif(1)
    e = R[Mar2[1,i], Mar2[2,i]]                                       #On créer des variables pour simplifier l'écriture
    f = R[Mar2[1,i], Mar2[2,i] + K]
    g = S[Mar2[1,i], Mar2[2,i] + K]
    h = S[Mar2[1,i], Mar2[2,i]]
    
    if((a > 10**300) | (b > 10**300) | (c > 10**300) | (d > 10**300)){
      warning('renforcement trop important')
      return(list(Mar1 = matrix(Mar1[Mar1 > 0], nrow = 2, byrow = FALSE), Mar2 = matrix(Mar2[Mar2 > 0], nrow = 2, byrow = FALSE), P = P, Q = Q, R = R, S = S, i = i))
    }
    
    if((e > 10**300) | (f > 10**300) | (g > 10**300) | (h > 10**300)){
      warning('renforcement trop important')
      return(list(Mar1 = matrix(Mar1[Mar1 > 0], nrow = 2, byrow = FALSE), Mar2 = matrix(Mar2[Mar2 > 0], nrow = 2, byrow = FALSE), P = P, Q = Q, R = R, S = S, i = i))
    }
    
    if(W < 0.5){
      if(U < (a / (a + b + c + d))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
        P[Mar1[1,i], Mar1[2,i]] = beta * P[Mar1[1,i], Mar1[2,i]] + alpha                           #On va a gauche et on ajoute delta 
        R[Mar1[1,i], Mar1[2,i]] = max(0, delta * R[Mar1[1,i], Mar1[2,i]] - gamma)
        Mar1[1, i + 1] = Mar1[1, i] - 1
        Mar1[2, i + 1] = Mar1[2, i]
      }
      else if(U < ((a + b) / (a + b + c + d))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
        P[Mar1[1,i], Mar1[2,i] + K] = beta * P[Mar1[1,i], Mar1[2,i] + K] + alpha                                               #On va a droite et on enleve delta
        R[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * R[Mar1[1,i], Mar1[2,i] + K] - gamma)
        Mar1[1, i + 1] = Mar1[1, i] + 1
        Mar1[2, i + 1] = Mar1[2, i]
      }
      else if(U < ((a + b + c) / (a + b + c + d))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
        Q[Mar1[1,i], Mar1[2,i] + K] = beta * Q[Mar1[1,i], Mar1[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
        S[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * S[Mar1[1,i], Mar1[2,i] + K] - gamma)
        Mar1[1, i + 1] = Mar1[1, i]
        Mar1[2, i + 1] = Mar1[2, i] + 1
      }
      else if((U != 1) | (Mar1[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
        Q[Mar1[1,i], Mar1[2,i]] = beta * Q[Mar1[1,i], Mar1[2,i]] + alpha                                                      #On va derriere et on enleve delta
        S[Mar1[1,i], Mar1[2,i]] = max(0, delta * S[Mar1[1,i], Mar1[2,i]] - gamma)
        Mar1[1, i + 1] = Mar1[1, i]
        Mar1[2, i + 1] = Mar1[2, i] - 1
      }
      else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
        Mar1[1, i + 1] = Mar1[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
        Mar1[2, i + 1] = Mar1[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
      } 
      
      if(V < (e / (e + f + g + h))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
        R[Mar2[1,i], Mar2[2,i]] = beta * R[Mar2[1,i], Mar2[2,i]] + alpha                           #On va a gauche et on ajoute delta 
        P[Mar2[1,i], Mar2[2,i]] = max(0, delta * P[Mar2[1,i], Mar2[2,i]] - gamma)
        Mar2[1, i + 1] = Mar2[1, i] - 1
        Mar2[2, i + 1] = Mar2[2, i]
      }
      else if(V < ((e + f) / (e + f + g + h))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
        R[Mar2[1,i], Mar2[2,i] + K] = beta * R[Mar2[1,i], Mar2[2,i] + K] + alpha                                               #On va a droite et on enleve delta
        P[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * P[Mar2[1,i],  Mar2[2,i] + K] - gamma)
        Mar2[1, i + 1] = Mar2[1, i] + 1
        Mar2[2, i + 1] = Mar2[2, i]
      }
      else if(V < ((e + f + g) / (e + f + g + h))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
        S[Mar2[1,i], Mar2[2,i] + K] = beta * S[Mar2[1,i], Mar2[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
        Q[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * Q[Mar2[1,i], Mar2[2,i] + K] - gamma)
        Mar2[1, i + 1] = Mar2[1, i]
        Mar2[2, i + 1] = Mar2[2, i] + 1
      }
      else if((V != 1) | (Mar2[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
        S[Mar2[1,i], Mar2[2,i]] = beta * S[Mar2[1,i], Mar2[2,i]] + alpha                                                      #On va derriere et on enleve delta
        Q[Mar2[1,i], Mar2[2,i]] = max(0, delta * Q[Mar2[1,i], Mar2[2,i]] - gamma)
        Mar2[1, i + 1] = Mar2[1, i]
        Mar2[2, i + 1] = Mar2[2, i] - 1
      }
      else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
        Mar2[1, i + 1] = Mar2[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
        Mar2[2, i + 1] = Mar2[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
      } 
    } else{
        if(V < (e / (e + f + g + h))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
          R[Mar2[1,i], Mar2[2,i]] = beta * R[Mar2[1,i], Mar2[2,i]] + alpha                           #On va a gauche et on ajoute delta 
          P[Mar2[1,i], Mar2[2,i]] = max(0, delta * P[Mar2[1,i], Mar2[2,i]] - gamma)
          Mar2[1, i + 1] = Mar2[1, i] - 1
          Mar2[2, i + 1] = Mar2[2, i]
        }
        else if(V < ((e + f) / (e + f + g + h))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
          R[Mar2[1,i], Mar2[2,i] + K] = beta * R[Mar2[1,i], Mar2[2,i] + K] + alpha                                               #On va a droite et on enleve delta
          P[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * P[Mar2[1,i],  Mar2[2,i] + K] - gamma)
          Mar2[1, i + 1] = Mar2[1, i] + 1
          Mar2[2, i + 1] = Mar2[2, i]
        }
        else if(V < ((e + f + g) / (e + f + g + h))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
          S[Mar2[1,i], Mar2[2,i] + K] = beta * S[Mar2[1,i], Mar2[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
          Q[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * Q[Mar2[1,i], Mar2[2,i] + K] - gamma)
          Mar2[1, i + 1] = Mar2[1, i]
          Mar2[2, i + 1] = Mar2[2, i] + 1
        }
        else if((V != 1) | (Mar2[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
          S[Mar2[1,i], Mar2[2,i]] = beta * S[Mar2[1,i], Mar2[2,i]] + alpha                                                      #On va derriere et on enleve delta
          Q[Mar2[1,i], Mar2[2,i]] = max(0, delta * Q[Mar2[1,i], Mar2[2,i]] - gamma)
          Mar2[1, i + 1] = Mar2[1, i]
          Mar2[2, i + 1] = Mar2[2, i] - 1
        }
        else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
          Mar2[1, i + 1] = Mar2[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
          Mar2[2, i + 1] = Mar2[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
        } 
        
        if(U < (a / (a + b + c + d))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
          P[Mar1[1,i], Mar1[2,i]] = beta * P[Mar1[1,i], Mar1[2,i]] + alpha                           #On va a gauche et on ajoute delta 
          R[Mar1[1,i], Mar1[2,i]] = max(0, delta * R[Mar1[1,i], Mar1[2,i]] - gamma)
          Mar1[1, i + 1] = Mar1[1, i] - 1
          Mar1[2, i + 1] = Mar1[2, i]
        }
        else if(U < ((a + b) / (a + b + c + d))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
          P[Mar1[1,i], Mar1[2,i] + K] = beta * P[Mar1[1,i], Mar1[2,i] + K] + alpha                                               #On va a droite et on enleve delta
          R[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * R[Mar1[1,i], Mar1[2,i] + K] - gamma)
          Mar1[1, i + 1] = Mar1[1, i] + 1
          Mar1[2, i + 1] = Mar1[2, i]
        }
        else if(U < ((a + b + c) / (a + b + c + d))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
          Q[Mar1[1,i], Mar1[2,i] + K] = beta * Q[Mar1[1,i], Mar1[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
          S[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * S[Mar1[1,i], Mar1[2,i] + K] - gamma)
          Mar1[1, i + 1] = Mar1[1, i]
          Mar1[2, i + 1] = Mar1[2, i] + 1
        }
        else if((U != 1) | (Mar1[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
          Q[Mar1[1,i], Mar1[2,i]] = beta * Q[Mar1[1,i], Mar1[2,i]] + alpha                                                      #On va derriere et on enleve delta
          S[Mar1[1,i], Mar1[2,i]] = max(0, delta * S[Mar1[1,i], Mar1[2,i]] - gamma)
          Mar1[1, i + 1] = Mar1[1, i]
          Mar1[2, i + 1] = Mar1[2, i] - 1
        }
        else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
          Mar1[1, i + 1] = Mar1[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
          Mar1[2, i + 1] = Mar1[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
        } 
      }
    }
  return(list(Mar1 = Mar1, Mar2 = Mar2, P = P,Q = Q,  R = R, S = S))
}


MAR3pop = function(P, Q, R, S, U, V, alpha, beta, gamma, delta, N, xdebut, ydebut, xdebut1, ydebut1, xdebut2, ydebut2){
  Mar1 = matrix(0, nrow = 2, ncol = N)
  Mar2 = matrix(0, nrow = 2, ncol = N)
  Mar3 = matrix(0, nrow = 2, ncol = N)
  Mar1[1,1] = xdebut
  Mar2[1,1] = xdebut1
  Mar3[1,1] = xdebut2
  Mar1[2,1] = ydebut
  Mar2[2,1] = ydebut1
  Mar3[2,1] = ydebut2
  K = nrow(P)
  for(i in 1:(N-1)){
    x = runif(1)
    
    u = runif(1)                                                  #On tire un nombre aleatoire entre 0 et 1
    a = P[Mar1[1,i], Mar1[2,i]]                                       #On créer des variables pour simplifier l'écriture
    b = P[Mar1[1,i], Mar1[2,i] + K]
    c = Q[Mar1[1,i], Mar1[2,i] + K]
    d = Q[Mar1[1,i], Mar1[2,i]]
    
    v = runif(1)
    e = R[Mar2[1,i], Mar2[2,i]]                                       #On créer des variables pour simplifier l'écriture
    f = R[Mar2[1,i], Mar2[2,i] + K]
    g = S[Mar2[1,i], Mar2[2,i] + K]
    h = S[Mar2[1,i], Mar2[2,i]]
    
    w = runif(1)
    j = U[Mar3[1,i], Mar3[2,i]]                                       #On créer des variables pour simplifier l'écriture
    k = U[Mar3[1,i], Mar3[2,i] + K]
    l = V[Mar3[1,i], Mar3[2,i] + K]
    m = V[Mar3[1,i], Mar3[2,i]]
    
    if((a > 10**300) | (b > 10**300) | (c > 10**300) | (d > 10**300)){
      warning('renforcement trop important')
      return(list(Mar1 = matrix(Mar1[Mar1 > 0], nrow = 2, byrow = FALSE), Mar2 = matrix(Mar2[Mar2 > 0], nrow = 2, byrow = FALSE), Mar3 = matrix(Mar3[Mar3 > 0], nrow = 2, byrow = FALSE), P = P, Q = Q, R = R, S = S, U = U, V = V, i = i))
    }
    
    if((e > 10**300) | (f > 10**300) | (g > 10**300) | (h > 10**300)){
      warning('renforcement trop important')
      return(list(Mar1 = matrix(Mar1[Mar1 > 0], nrow = 2, byrow = FALSE), Mar2 = matrix(Mar2[Mar2 > 0], nrow = 2, byrow = FALSE), Mar3 = matrix(Mar3[Mar3 > 0], nrow = 2, byrow = FALSE), P = P, Q = Q, R = R, S = S, U = U, V = V, i = i))
    }
    
    if((j > 10**300) | (k > 10**300) | (l > 10**300) | (m > 10**300)){
      warning('renforcement trop important')
      return(list(Mar1 = matrix(Mar1[Mar1 > 0], nrow = 2, byrow = FALSE), Mar2 = matrix(Mar2[Mar2 > 0], nrow = 2, byrow = FALSE), Mar3 = matrix(Mar3[Mar3 > 0], nrow = 2, byrow = FALSE), P = P, Q = Q, R = R, S = S, U = U, V = V, i = i))
    }
    
    if(x < (1/3)){
      
      if(u < (a / (a + b + c + d))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
        P[Mar1[1,i], Mar1[2,i]] = beta * P[Mar1[1,i], Mar1[2,i]] + alpha                           #On va a gauche et on ajoute delta 
        R[Mar1[1,i], Mar1[2,i]] = max(0, delta * R[Mar1[1,i], Mar1[2,i]] - gamma)
        U[Mar1[1,i], Mar1[2,i]] = max(0, delta * U[Mar1[1,i], Mar1[2,i]] - gamma)
        Mar1[1, i + 1] = Mar1[1, i] - 1
        Mar1[2, i + 1] = Mar1[2, i]
      }
      else if(u < ((a + b) / (a + b + c + d))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
        P[Mar1[1,i], Mar1[2,i] + K] = beta * P[Mar1[1,i], Mar1[2,i] + K] + alpha                                               #On va a droite et on enleve delta
        R[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * R[Mar1[1,i], Mar1[2,i] + K] - gamma)
        U[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * U[Mar1[1,i], Mar1[2,i] + K] - gamma)
        Mar1[1, i + 1] = Mar1[1, i] + 1
        Mar1[2, i + 1] = Mar1[2, i]
      }
      else if(u < ((a + b + c) / (a + b + c + d))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
        Q[Mar1[1,i], Mar1[2,i] + K] = beta * Q[Mar1[1,i], Mar1[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
        S[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * S[Mar1[1,i], Mar1[2,i] + K] - gamma)
        V[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * V[Mar1[1,i], Mar1[2,i] + K] - gamma)
        Mar1[1, i + 1] = Mar1[1, i]
        Mar1[2, i + 1] = Mar1[2, i] + 1
      }
      else if((u != 1) | (Mar1[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
        Q[Mar1[1,i], Mar1[2,i]] = beta * Q[Mar1[1,i], Mar1[2,i]] + alpha                                                      #On va derriere et on enleve delta
        S[Mar1[1,i], Mar1[2,i]] = max(0, delta * S[Mar1[1,i], Mar1[2,i]] - gamma)
        V[Mar1[1,i], Mar1[2,i]] = max(0, delta * V[Mar1[1,i], Mar1[2,i]] - gamma)
        Mar1[1, i + 1] = Mar1[1, i]
        Mar1[2, i + 1] = Mar1[2, i] - 1
      }
      else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
        Mar1[1, i + 1] = Mar1[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
        Mar1[2, i + 1] = Mar1[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
      }  
      if(x < 1/6){
      
        if(v < (e / (e + f + g + h))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
          R[Mar2[1,i], Mar2[2,i]] = beta * R[Mar2[1,i], Mar2[2,i]] + alpha                           #On va a gauche et on ajoute delta 
          P[Mar2[1,i], Mar2[2,i]] = max(0, delta * P[Mar2[1,i], Mar2[2,i]] - gamma)
          U[Mar2[1,i], Mar2[2,i]] = max(0, delta * U[Mar2[1,i], Mar2[2,i]] - gamma)
          Mar2[1, i + 1] = Mar2[1, i] - 1
          Mar2[2, i + 1] = Mar2[2, i]
        }
        else if(v < ((e + f) / (e + f + g + h))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
          R[Mar2[1,i], Mar2[2,i] + K] = beta * R[Mar2[1,i], Mar2[2,i] + K] + alpha                                               #On va a droite et on enleve delta
          P[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * P[Mar2[1,i],  Mar2[2,i] + K] - gamma)
          U[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * U[Mar2[1,i],  Mar2[2,i] + K] - gamma)
          Mar2[1, i + 1] = Mar2[1, i] + 1
          Mar2[2, i + 1] = Mar2[2, i]
        }
        else if(v < ((e + f + g) / (e + f + g + h))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
          S[Mar2[1,i], Mar2[2,i] + K] = beta * S[Mar2[1,i], Mar2[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
          Q[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * Q[Mar2[1,i], Mar2[2,i] + K] - gamma)
          V[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * V[Mar2[1,i], Mar2[2,i] + K] - gamma)
          Mar2[1, i + 1] = Mar2[1, i]
          Mar2[2, i + 1] = Mar2[2, i] + 1
        }
        else if((v != 1) | (Mar2[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
          S[Mar2[1,i], Mar2[2,i]] = beta * S[Mar2[1,i], Mar2[2,i]] + alpha                                                      #On va derriere et on enleve delta
          Q[Mar2[1,i], Mar2[2,i]] = max(0, delta * Q[Mar2[1,i], Mar2[2,i]] - gamma)
          V[Mar2[1,i], Mar2[2,i]] = max(0, delta * V[Mar2[1,i], Mar2[2,i]] - gamma)
          Mar2[1, i + 1] = Mar2[1, i]
          Mar2[2, i + 1] = Mar2[2, i] - 1
        }
        else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
          Mar2[1, i + 1] = Mar2[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
          Mar2[2, i + 1] = Mar2[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
        }
        
        if(w < (j / (j + k + l + m))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
          U[Mar3[1,i], Mar3[2,i]] = beta * U[Mar3[1,i], Mar3[2,i]] + alpha                           #On va a gauche et on ajoute delta 
          P[Mar3[1,i], Mar3[2,i]] = max(0, delta * P[Mar3[1,i], Mar3[2,i]] - gamma)
          R[Mar3[1,i], Mar3[2,i]] = max(0, delta * R[Mar3[1,i], Mar3[2,i]] - gamma)
          Mar3[1, i + 1] = Mar3[1, i] - 1
          Mar3[2, i + 1] = Mar3[2, i]
        }
        else if(w < ((j + k) / (j + k + l + m))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
          U[Mar3[1,i], Mar3[2,i] + K] = beta * U[Mar3[1,i], Mar3[2,i] + K] + alpha                                               #On va a droite et on enleve delta
          P[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * P[Mar3[1,i],  Mar3[2,i] + K] - gamma)
          R[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * R[Mar3[1,i],  Mar3[2,i] + K] - gamma)
          Mar3[1, i + 1] = Mar3[1, i] + 1
          Mar3[2, i + 1] = Mar3[2, i]
        }
        else if(w < ((j + k + l) / (j + k + l + m))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
          V[Mar3[1,i], Mar3[2,i] + K] = beta * V[Mar3[1,i], Mar3[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
          Q[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * Q[Mar3[1,i], Mar3[2,i] + K] - gamma)
          S[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * S[Mar3[1,i], Mar3[2,i] + K] - gamma)
          Mar3[1, i + 1] = Mar3[1, i]
          Mar3[2, i + 1] = Mar3[2, i] + 1
        }
        else if((w != 1) | (Mar3[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
          V[Mar3[1,i], Mar3[2,i]] = beta * V[Mar3[1,i], Mar3[2,i]] + alpha                                                      #On va derriere et on enleve delta
          Q[Mar3[1,i], Mar3[2,i]] = max(0, delta * Q[Mar3[1,i], Mar3[2,i]] - gamma)
          S[Mar3[1,i], Mar3[2,i]] = max(0, delta * S[Mar3[1,i], Mar3[2,i]] - gamma)
          Mar3[1, i + 1] = Mar3[1, i]
          Mar3[2, i + 1] = Mar3[2, i] - 1
        }
        else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
          Mar3[1, i + 1] = Mar3[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
          Mar3[2, i + 1] = Mar3[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
        }
      }
      else{
        
        if(w < (j / (j + k + l + m))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
          U[Mar3[1,i], Mar3[2,i]] = beta * U[Mar3[1,i], Mar3[2,i]] + alpha                           #On va a gauche et on ajoute delta 
          P[Mar3[1,i], Mar3[2,i]] = max(0, delta * P[Mar3[1,i], Mar3[2,i]] - gamma)
          R[Mar3[1,i], Mar3[2,i]] = max(0, delta * R[Mar3[1,i], Mar3[2,i]] - gamma)
          Mar3[1, i + 1] = Mar3[1, i] - 1
          Mar3[2, i + 1] = Mar3[2, i]
        }
        else if(w < ((j + k) / (j + k + l + m))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
          U[Mar3[1,i], Mar3[2,i] + K] = beta * U[Mar3[1,i], Mar3[2,i] + K] + alpha                                               #On va a droite et on enleve delta
          P[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * P[Mar3[1,i],  Mar3[2,i] + K] - gamma)
          R[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * R[Mar3[1,i],  Mar3[2,i] + K] - gamma)
          Mar3[1, i + 1] = Mar3[1, i] + 1
          Mar3[2, i + 1] = Mar3[2, i]
        }
        else if(w < ((j + k + l) / (j + k + l + m))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
          V[Mar3[1,i], Mar3[2,i] + K] = beta * V[Mar3[1,i], Mar3[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
          Q[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * Q[Mar3[1,i], Mar3[2,i] + K] - gamma)
          S[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * S[Mar3[1,i], Mar3[2,i] + K] - gamma)
          Mar3[1, i + 1] = Mar3[1, i]
          Mar3[2, i + 1] = Mar3[2, i] + 1
        }
        else if((w != 1) | (Mar3[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
          V[Mar3[1,i], Mar3[2,i]] = beta * V[Mar3[1,i], Mar3[2,i]] + alpha                                                      #On va derriere et on enleve delta
          Q[Mar3[1,i], Mar3[2,i]] = max(0, delta * Q[Mar3[1,i], Mar3[2,i]] - gamma)
          S[Mar3[1,i], Mar3[2,i]] = max(0, delta * S[Mar3[1,i], Mar3[2,i]] - gamma)
          Mar3[1, i + 1] = Mar3[1, i]
          Mar3[2, i + 1] = Mar3[2, i] - 1
        }
        else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
          Mar3[1, i + 1] = Mar3[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
          Mar3[2, i + 1] = Mar3[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
        }
        
        if(v < (e / (e + f + g + h))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
          R[Mar2[1,i], Mar2[2,i]] = beta * R[Mar2[1,i], Mar2[2,i]] + alpha                           #On va a gauche et on ajoute delta 
          P[Mar2[1,i], Mar2[2,i]] = max(0, delta * P[Mar2[1,i], Mar2[2,i]] - gamma)
          U[Mar2[1,i], Mar2[2,i]] = max(0, delta * U[Mar2[1,i], Mar2[2,i]] - gamma)
          Mar2[1, i + 1] = Mar2[1, i] - 1
          Mar2[2, i + 1] = Mar2[2, i]
        }
        else if(v < ((e + f) / (e + f + g + h))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
          R[Mar2[1,i], Mar2[2,i] + K] = beta * R[Mar2[1,i], Mar2[2,i] + K] + alpha                                               #On va a droite et on enleve delta
          P[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * P[Mar2[1,i],  Mar2[2,i] + K] - gamma)
          U[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * U[Mar2[1,i],  Mar2[2,i] + K] - gamma)
          Mar2[1, i + 1] = Mar2[1, i] + 1
          Mar2[2, i + 1] = Mar2[2, i]
        }
        else if(v < ((e + f + g) / (e + f + g + h))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
          S[Mar2[1,i], Mar2[2,i] + K] = beta * S[Mar2[1,i], Mar2[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
          Q[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * Q[Mar2[1,i], Mar2[2,i] + K] - gamma)
          V[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * V[Mar2[1,i], Mar2[2,i] + K] - gamma)
          Mar2[1, i + 1] = Mar2[1, i]
          Mar2[2, i + 1] = Mar2[2, i] + 1
        }
        else if((v != 1) | (Mar2[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
          S[Mar2[1,i], Mar2[2,i]] = beta * S[Mar2[1,i], Mar2[2,i]] + alpha                                                      #On va derriere et on enleve delta
          Q[Mar2[1,i], Mar2[2,i]] = max(0, delta * Q[Mar2[1,i], Mar2[2,i]] - gamma)
          V[Mar2[1,i], Mar2[2,i]] = max(0, delta * V[Mar2[1,i], Mar2[2,i]] - gamma)
          Mar2[1, i + 1] = Mar2[1, i]
          Mar2[2, i + 1] = Mar2[2, i] - 1
        }
        else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
          Mar2[1, i + 1] = Mar2[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
          Mar2[2, i + 1] = Mar2[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
        }
        
      }
    }
    else if(x < (2/3)){
        
      if(v < (e / (e + f + g + h))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
        R[Mar2[1,i], Mar2[2,i]] = beta * R[Mar2[1,i], Mar2[2,i]] + alpha                           #On va a gauche et on ajoute delta 
        P[Mar2[1,i], Mar2[2,i]] = max(0, delta * P[Mar2[1,i], Mar2[2,i]] - gamma)
        U[Mar2[1,i], Mar2[2,i]] = max(0, delta * U[Mar2[1,i], Mar2[2,i]] - gamma)
        Mar2[1, i + 1] = Mar2[1, i] - 1
        Mar2[2, i + 1] = Mar2[2, i]
      }
      else if(v < ((e + f) / (e + f + g + h))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
        R[Mar2[1,i], Mar2[2,i] + K] = beta * R[Mar2[1,i], Mar2[2,i] + K] + alpha                                               #On va a droite et on enleve delta
        P[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * P[Mar2[1,i],  Mar2[2,i] + K] - gamma)
        U[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * U[Mar2[1,i],  Mar2[2,i] + K] - gamma)
        Mar2[1, i + 1] = Mar2[1, i] + 1
        Mar2[2, i + 1] = Mar2[2, i]
      }
      else if(v < ((e + f + g) / (e + f + g + h))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
        S[Mar2[1,i], Mar2[2,i] + K] = beta * S[Mar2[1,i], Mar2[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
        Q[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * Q[Mar2[1,i], Mar2[2,i] + K] - gamma)
        V[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * V[Mar2[1,i], Mar2[2,i] + K] - gamma)
        Mar2[1, i + 1] = Mar2[1, i]
        Mar2[2, i + 1] = Mar2[2, i] + 1
      }
      else if((v != 1) | (Mar2[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
        S[Mar2[1,i], Mar2[2,i]] = beta * S[Mar2[1,i], Mar2[2,i]] + alpha                                                      #On va derriere et on enleve delta
        Q[Mar2[1,i], Mar2[2,i]] = max(0, delta * Q[Mar2[1,i], Mar2[2,i]] - gamma)
        V[Mar2[1,i], Mar2[2,i]] = max(0, delta * V[Mar2[1,i], Mar2[2,i]] - gamma)
        Mar2[1, i + 1] = Mar2[1, i]
        Mar2[2, i + 1] = Mar2[2, i] - 1
      }
      else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
        Mar2[1, i + 1] = Mar2[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
        Mar2[2, i + 1] = Mar2[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
      }
      
      if(x < 1/2){
        
        if(u < (a / (a + b + c + d))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
          P[Mar1[1,i], Mar1[2,i]] = beta * P[Mar1[1,i], Mar1[2,i]] + alpha                           #On va a gauche et on ajoute delta 
          R[Mar1[1,i], Mar1[2,i]] = max(0, delta * R[Mar1[1,i], Mar1[2,i]] - gamma)
          U[Mar1[1,i], Mar1[2,i]] = max(0, delta * U[Mar1[1,i], Mar1[2,i]] - gamma)
          Mar1[1, i + 1] = Mar1[1, i] - 1
          Mar1[2, i + 1] = Mar1[2, i]
        }
        else if(u < ((a + b) / (a + b + c + d))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
          P[Mar1[1,i], Mar1[2,i] + K] = beta * P[Mar1[1,i], Mar1[2,i] + K] + alpha                                               #On va a droite et on enleve delta
          R[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * R[Mar1[1,i], Mar1[2,i] + K] - gamma)
          U[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * U[Mar1[1,i], Mar1[2,i] + K] - gamma)
          Mar1[1, i + 1] = Mar1[1, i] + 1
          Mar1[2, i + 1] = Mar1[2, i]
        }
        else if(u < ((a + b + c) / (a + b + c + d))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
          Q[Mar1[1,i], Mar1[2,i] + K] = beta * Q[Mar1[1,i], Mar1[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
          S[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * S[Mar1[1,i], Mar1[2,i] + K] - gamma)
          V[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * V[Mar1[1,i], Mar1[2,i] + K] - gamma)
          Mar1[1, i + 1] = Mar1[1, i]
          Mar1[2, i + 1] = Mar1[2, i] + 1
        }
        else if((u != 1) | (Mar1[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
          Q[Mar1[1,i], Mar1[2,i]] = beta * Q[Mar1[1,i], Mar1[2,i]] + alpha                                                      #On va derriere et on enleve delta
          S[Mar1[1,i], Mar1[2,i]] = max(0, delta * S[Mar1[1,i], Mar1[2,i]] - gamma)
          V[Mar1[1,i], Mar1[2,i]] = max(0, delta * V[Mar1[1,i], Mar1[2,i]] - gamma)
          Mar1[1, i + 1] = Mar1[1, i]
          Mar1[2, i + 1] = Mar1[2, i] - 1
        }
        else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
          Mar1[1, i + 1] = Mar1[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
          Mar1[2, i + 1] = Mar1[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
        }  
        
        if(w < (j / (j + k + l + m))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
          U[Mar3[1,i], Mar3[2,i]] = beta * U[Mar3[1,i], Mar3[2,i]] + alpha                           #On va a gauche et on ajoute delta 
          P[Mar3[1,i], Mar3[2,i]] = max(0, delta * P[Mar3[1,i], Mar3[2,i]] - gamma)
          R[Mar3[1,i], Mar3[2,i]] = max(0, delta * R[Mar3[1,i], Mar3[2,i]] - gamma)
          Mar3[1, i + 1] = Mar3[1, i] - 1
          Mar3[2, i + 1] = Mar3[2, i]
        }
        else if(w < ((j + k) / (j + k + l + m))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
          U[Mar3[1,i], Mar3[2,i] + K] = beta * U[Mar3[1,i], Mar3[2,i] + K] + alpha                                               #On va a droite et on enleve delta
          P[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * P[Mar3[1,i],  Mar3[2,i] + K] - gamma)
          R[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * R[Mar3[1,i],  Mar3[2,i] + K] - gamma)
          Mar3[1, i + 1] = Mar3[1, i] + 1
          Mar3[2, i + 1] = Mar3[2, i]
        }
        else if(w < ((j + k + l) / (j + k + l + m))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
          V[Mar3[1,i], Mar3[2,i] + K] = beta * V[Mar3[1,i], Mar3[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
          Q[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * Q[Mar3[1,i], Mar3[2,i] + K] - gamma)
          S[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * S[Mar3[1,i], Mar3[2,i] + K] - gamma)
          Mar3[1, i + 1] = Mar3[1, i]
          Mar3[2, i + 1] = Mar3[2, i] + 1
        }
        else if((w != 1) | (Mar3[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
          V[Mar3[1,i], Mar3[2,i]] = beta * V[Mar3[1,i], Mar3[2,i]] + alpha                                                      #On va derriere et on enleve delta
          Q[Mar3[1,i], Mar3[2,i]] = max(0, delta * Q[Mar3[1,i], Mar3[2,i]] - gamma)
          S[Mar3[1,i], Mar3[2,i]] = max(0, delta * S[Mar3[1,i], Mar3[2,i]] - gamma)
          Mar3[1, i + 1] = Mar3[1, i]
          Mar3[2, i + 1] = Mar3[2, i] - 1
        }
        else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
          Mar3[1, i + 1] = Mar3[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
          Mar3[2, i + 1] = Mar3[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
        }
      }
      
      else{
        
        if(w < (j / (j + k + l + m))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
          U[Mar3[1,i], Mar3[2,i]] = beta * U[Mar3[1,i], Mar3[2,i]] + alpha                           #On va a gauche et on ajoute delta 
          P[Mar3[1,i], Mar3[2,i]] = max(0, delta * P[Mar3[1,i], Mar3[2,i]] - gamma)
          R[Mar3[1,i], Mar3[2,i]] = max(0, delta * R[Mar3[1,i], Mar3[2,i]] - gamma)
          Mar3[1, i + 1] = Mar3[1, i] - 1
          Mar3[2, i + 1] = Mar3[2, i]
        }
        else if(w < ((j + k) / (j + k + l + m))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
          U[Mar3[1,i], Mar3[2,i] + K] = beta * U[Mar3[1,i], Mar3[2,i] + K] + alpha                                               #On va a droite et on enleve delta
          P[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * P[Mar3[1,i],  Mar3[2,i] + K] - gamma)
          R[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * R[Mar3[1,i],  Mar3[2,i] + K] - gamma)
          Mar3[1, i + 1] = Mar3[1, i] + 1
          Mar3[2, i + 1] = Mar3[2, i]
        }
        else if(w < ((j + k + l) / (j + k + l + m))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
          V[Mar3[1,i], Mar3[2,i] + K] = beta * V[Mar3[1,i], Mar3[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
          Q[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * Q[Mar3[1,i], Mar3[2,i] + K] - gamma)
          S[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * S[Mar3[1,i], Mar3[2,i] + K] - gamma)
          Mar3[1, i + 1] = Mar3[1, i]
          Mar3[2, i + 1] = Mar3[2, i] + 1
        }
        else if((w != 1) | (Mar3[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
          V[Mar3[1,i], Mar3[2,i]] = beta * V[Mar3[1,i], Mar3[2,i]] + alpha                                                      #On va derriere et on enleve delta
          Q[Mar3[1,i], Mar3[2,i]] = max(0, delta * Q[Mar3[1,i], Mar3[2,i]] - gamma)
          S[Mar3[1,i], Mar3[2,i]] = max(0, delta * S[Mar3[1,i], Mar3[2,i]] - gamma)
          Mar3[1, i + 1] = Mar3[1, i]
          Mar3[2, i + 1] = Mar3[2, i] - 1
        }
        else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
          Mar3[1, i + 1] = Mar3[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
          Mar3[2, i + 1] = Mar3[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
        }
        
        if(u < (a / (a + b + c + d))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
          P[Mar1[1,i], Mar1[2,i]] = beta * P[Mar1[1,i], Mar1[2,i]] + alpha                           #On va a gauche et on ajoute delta 
          R[Mar1[1,i], Mar1[2,i]] = max(0, delta * R[Mar1[1,i], Mar1[2,i]] - gamma)
          U[Mar1[1,i], Mar1[2,i]] = max(0, delta * U[Mar1[1,i], Mar1[2,i]] - gamma)
          Mar1[1, i + 1] = Mar1[1, i] - 1
          Mar1[2, i + 1] = Mar1[2, i]
        }
        else if(u < ((a + b) / (a + b + c + d))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
          P[Mar1[1,i], Mar1[2,i] + K] = beta * P[Mar1[1,i], Mar1[2,i] + K] + alpha                                               #On va a droite et on enleve delta
          R[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * R[Mar1[1,i], Mar1[2,i] + K] - gamma)
          U[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * U[Mar1[1,i], Mar1[2,i] + K] - gamma)
          Mar1[1, i + 1] = Mar1[1, i] + 1
          Mar1[2, i + 1] = Mar1[2, i]
        }
        else if(u < ((a + b + c) / (a + b + c + d))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
          Q[Mar1[1,i], Mar1[2,i] + K] = beta * Q[Mar1[1,i], Mar1[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
          S[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * S[Mar1[1,i], Mar1[2,i] + K] - gamma)
          V[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * V[Mar1[1,i], Mar1[2,i] + K] - gamma)
          Mar1[1, i + 1] = Mar1[1, i]
          Mar1[2, i + 1] = Mar1[2, i] + 1
        }
        else if((u != 1) | (Mar1[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
          Q[Mar1[1,i], Mar1[2,i]] = beta * Q[Mar1[1,i], Mar1[2,i]] + alpha                                                      #On va derriere et on enleve delta
          S[Mar1[1,i], Mar1[2,i]] = max(0, delta * S[Mar1[1,i], Mar1[2,i]] - gamma)
          V[Mar1[1,i], Mar1[2,i]] = max(0, delta * V[Mar1[1,i], Mar1[2,i]] - gamma)
          Mar1[1, i + 1] = Mar1[1, i]
          Mar1[2, i + 1] = Mar1[2, i] - 1
        }
        else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
          Mar1[1, i + 1] = Mar1[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
          Mar1[2, i + 1] = Mar1[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
        }  
        
        
      }
      
    }
    else{
      
      if(w < (j / (j + k + l + m))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
        U[Mar3[1,i], Mar3[2,i]] = beta * U[Mar3[1,i], Mar3[2,i]] + alpha                           #On va a gauche et on ajoute delta 
        P[Mar3[1,i], Mar3[2,i]] = max(0, delta * P[Mar3[1,i], Mar3[2,i]] - gamma)
        R[Mar3[1,i], Mar3[2,i]] = max(0, delta * R[Mar3[1,i], Mar3[2,i]] - gamma)
        Mar3[1, i + 1] = Mar3[1, i] - 1
        Mar3[2, i + 1] = Mar3[2, i]
      }
      else if(w < ((j + k) / (j + k + l + m))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
        U[Mar3[1,i], Mar3[2,i] + K] = beta * U[Mar3[1,i], Mar3[2,i] + K] + alpha                                               #On va a droite et on enleve delta
        P[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * P[Mar3[1,i],  Mar3[2,i] + K] - gamma)
        R[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * R[Mar3[1,i],  Mar3[2,i] + K] - gamma)
        Mar3[1, i + 1] = Mar3[1, i] + 1
        Mar3[2, i + 1] = Mar3[2, i]
      }
      else if(w < ((j + k + l) / (j + k + l + m))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
        V[Mar3[1,i], Mar3[2,i] + K] = beta * V[Mar3[1,i], Mar3[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
        Q[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * Q[Mar3[1,i], Mar3[2,i] + K] - gamma)
        S[Mar3[1,i], Mar3[2,i] + K] = max(0, delta * S[Mar3[1,i], Mar3[2,i] + K] - gamma)
        Mar3[1, i + 1] = Mar3[1, i]
        Mar3[2, i + 1] = Mar3[2, i] + 1
      }
      else if((w != 1) | (Mar3[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
        V[Mar3[1,i], Mar3[2,i]] = beta * V[Mar3[1,i], Mar3[2,i]] + alpha                                                      #On va derriere et on enleve delta
        Q[Mar3[1,i], Mar3[2,i]] = max(0, delta * Q[Mar3[1,i], Mar3[2,i]] - gamma)
        S[Mar3[1,i], Mar3[2,i]] = max(0, delta * S[Mar3[1,i], Mar3[2,i]] - gamma)
        Mar3[1, i + 1] = Mar3[1, i]
        Mar3[2, i + 1] = Mar3[2, i] - 1
      }
      else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
        Mar3[1, i + 1] = Mar3[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
        Mar3[2, i + 1] = Mar3[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
      }
      
      if(x < 5/6){
        
        if(u < (a / (a + b + c + d))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
          P[Mar1[1,i], Mar1[2,i]] = beta * P[Mar1[1,i], Mar1[2,i]] + alpha                           #On va a gauche et on ajoute delta 
          R[Mar1[1,i], Mar1[2,i]] = max(0, delta * R[Mar1[1,i], Mar1[2,i]] - gamma)
          U[Mar1[1,i], Mar1[2,i]] = max(0, delta * U[Mar1[1,i], Mar1[2,i]] - gamma)
          Mar1[1, i + 1] = Mar1[1, i] - 1
          Mar1[2, i + 1] = Mar1[2, i]
        }
        else if(u < ((a + b) / (a + b + c + d))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
          P[Mar1[1,i], Mar1[2,i] + K] = beta * P[Mar1[1,i], Mar1[2,i] + K] + alpha                                               #On va a droite et on enleve delta
          R[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * R[Mar1[1,i], Mar1[2,i] + K] - gamma)
          U[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * U[Mar1[1,i], Mar1[2,i] + K] - gamma)
          Mar1[1, i + 1] = Mar1[1, i] + 1
          Mar1[2, i + 1] = Mar1[2, i]
        }
        else if(u < ((a + b + c) / (a + b + c + d))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
          Q[Mar1[1,i], Mar1[2,i] + K] = beta * Q[Mar1[1,i], Mar1[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
          S[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * S[Mar1[1,i], Mar1[2,i] + K] - gamma)
          V[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * V[Mar1[1,i], Mar1[2,i] + K] - gamma)
          Mar1[1, i + 1] = Mar1[1, i]
          Mar1[2, i + 1] = Mar1[2, i] + 1
        }
        else if((u != 1) | (Mar1[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
          Q[Mar1[1,i], Mar1[2,i]] = beta * Q[Mar1[1,i], Mar1[2,i]] + alpha                                                      #On va derriere et on enleve delta
          S[Mar1[1,i], Mar1[2,i]] = max(0, delta * S[Mar1[1,i], Mar1[2,i]] - gamma)
          V[Mar1[1,i], Mar1[2,i]] = max(0, delta * V[Mar1[1,i], Mar1[2,i]] - gamma)
          Mar1[1, i + 1] = Mar1[1, i]
          Mar1[2, i + 1] = Mar1[2, i] - 1
        }
        else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
          Mar1[1, i + 1] = Mar1[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
          Mar1[2, i + 1] = Mar1[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
        } 
        
        if(v < (e / (e + f + g + h))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
          R[Mar2[1,i], Mar2[2,i]] = beta * R[Mar2[1,i], Mar2[2,i]] + alpha                           #On va a gauche et on ajoute delta 
          P[Mar2[1,i], Mar2[2,i]] = max(0, delta * P[Mar2[1,i], Mar2[2,i]] - gamma)
          U[Mar2[1,i], Mar2[2,i]] = max(0, delta * U[Mar2[1,i], Mar2[2,i]] - gamma)
          Mar2[1, i + 1] = Mar2[1, i] - 1
          Mar2[2, i + 1] = Mar2[2, i]
        }
        else if(v < ((e + f) / (e + f + g + h))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
          R[Mar2[1,i], Mar2[2,i] + K] = beta * R[Mar2[1,i], Mar2[2,i] + K] + alpha                                               #On va a droite et on enleve delta
          P[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * P[Mar2[1,i],  Mar2[2,i] + K] - gamma)
          U[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * U[Mar2[1,i],  Mar2[2,i] + K] - gamma)
          Mar2[1, i + 1] = Mar2[1, i] + 1
          Mar2[2, i + 1] = Mar2[2, i]
        }
        else if(v < ((e + f + g) / (e + f + g + h))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
          S[Mar2[1,i], Mar2[2,i] + K] = beta * S[Mar2[1,i], Mar2[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
          Q[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * Q[Mar2[1,i], Mar2[2,i] + K] - gamma)
          V[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * V[Mar2[1,i], Mar2[2,i] + K] - gamma)
          Mar2[1, i + 1] = Mar2[1, i]
          Mar2[2, i + 1] = Mar2[2, i] + 1
        }
        else if((v != 1) | (Mar2[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
          S[Mar2[1,i], Mar2[2,i]] = beta * S[Mar2[1,i], Mar2[2,i]] + alpha                                                      #On va derriere et on enleve delta
          Q[Mar2[1,i], Mar2[2,i]] = max(0, delta * Q[Mar2[1,i], Mar2[2,i]] - gamma)
          V[Mar2[1,i], Mar2[2,i]] = max(0, delta * V[Mar2[1,i], Mar2[2,i]] - gamma)
          Mar2[1, i + 1] = Mar2[1, i]
          Mar2[2, i + 1] = Mar2[2, i] - 1
        }
        else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
          Mar2[1, i + 1] = Mar2[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
          Mar2[2, i + 1] = Mar2[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
        }
        
      }
      else{
        
        if(v < (e / (e + f + g + h))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
          R[Mar2[1,i], Mar2[2,i]] = beta * R[Mar2[1,i], Mar2[2,i]] + alpha                           #On va a gauche et on ajoute delta 
          P[Mar2[1,i], Mar2[2,i]] = max(0, delta * P[Mar2[1,i], Mar2[2,i]] - gamma)
          U[Mar2[1,i], Mar2[2,i]] = max(0, delta * U[Mar2[1,i], Mar2[2,i]] - gamma)
          Mar2[1, i + 1] = Mar2[1, i] - 1
          Mar2[2, i + 1] = Mar2[2, i]
        }
        else if(v < ((e + f) / (e + f + g + h))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
          R[Mar2[1,i], Mar2[2,i] + K] = beta * R[Mar2[1,i], Mar2[2,i] + K] + alpha                                               #On va a droite et on enleve delta
          P[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * P[Mar2[1,i],  Mar2[2,i] + K] - gamma)
          U[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * U[Mar2[1,i],  Mar2[2,i] + K] - gamma)
          Mar2[1, i + 1] = Mar2[1, i] + 1
          Mar2[2, i + 1] = Mar2[2, i]
        }
        else if(v < ((e + f + g) / (e + f + g + h))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
          S[Mar2[1,i], Mar2[2,i] + K] = beta * S[Mar2[1,i], Mar2[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
          Q[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * Q[Mar2[1,i], Mar2[2,i] + K] - gamma)
          V[Mar2[1,i], Mar2[2,i] + K] = max(0, delta * V[Mar2[1,i], Mar2[2,i] + K] - gamma)
          Mar2[1, i + 1] = Mar2[1, i]
          Mar2[2, i + 1] = Mar2[2, i] + 1
        }
        else if((v != 1) | (Mar2[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
          S[Mar2[1,i], Mar2[2,i]] = beta * S[Mar2[1,i], Mar2[2,i]] + alpha                                                      #On va derriere et on enleve delta
          Q[Mar2[1,i], Mar2[2,i]] = max(0, delta * Q[Mar2[1,i], Mar2[2,i]] - gamma)
          V[Mar2[1,i], Mar2[2,i]] = max(0, delta * V[Mar2[1,i], Mar2[2,i]] - gamma)
          Mar2[1, i + 1] = Mar2[1, i]
          Mar2[2, i + 1] = Mar2[2, i] - 1
        }
        else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
          Mar2[1, i + 1] = Mar2[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
          Mar2[2, i + 1] = Mar2[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
        }
        
        if(u < (a / (a + b + c + d))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
          P[Mar1[1,i], Mar1[2,i]] = beta * P[Mar1[1,i], Mar1[2,i]] + alpha                           #On va a gauche et on ajoute delta 
          R[Mar1[1,i], Mar1[2,i]] = max(0, delta * R[Mar1[1,i], Mar1[2,i]] - gamma)
          U[Mar1[1,i], Mar1[2,i]] = max(0, delta * U[Mar1[1,i], Mar1[2,i]] - gamma)
          Mar1[1, i + 1] = Mar1[1, i] - 1
          Mar1[2, i + 1] = Mar1[2, i]
        }
        else if(u < ((a + b) / (a + b + c + d))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
          P[Mar1[1,i], Mar1[2,i] + K] = beta * P[Mar1[1,i], Mar1[2,i] + K] + alpha                                               #On va a droite et on enleve delta
          R[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * R[Mar1[1,i], Mar1[2,i] + K] - gamma)
          U[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * U[Mar1[1,i], Mar1[2,i] + K] - gamma)
          Mar1[1, i + 1] = Mar1[1, i] + 1
          Mar1[2, i + 1] = Mar1[2, i]
        }
        else if(u < ((a + b + c) / (a + b + c + d))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
          Q[Mar1[1,i], Mar1[2,i] + K] = beta * Q[Mar1[1,i], Mar1[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
          S[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * S[Mar1[1,i], Mar1[2,i] + K] - gamma)
          V[Mar1[1,i], Mar1[2,i] + K] = max(0, delta * V[Mar1[1,i], Mar1[2,i] + K] - gamma)
          Mar1[1, i + 1] = Mar1[1, i]
          Mar1[2, i + 1] = Mar1[2, i] + 1
        }
        else if((u != 1) | (Mar1[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
          Q[Mar1[1,i], Mar1[2,i]] = beta * Q[Mar1[1,i], Mar1[2,i]] + alpha                                                      #On va derriere et on enleve delta
          S[Mar1[1,i], Mar1[2,i]] = max(0, delta * S[Mar1[1,i], Mar1[2,i]] - gamma)
          V[Mar1[1,i], Mar1[2,i]] = max(0, delta * V[Mar1[1,i], Mar1[2,i]] - gamma)
          Mar1[1, i + 1] = Mar1[1, i]
          Mar1[2, i + 1] = Mar1[2, i] - 1
        }
        else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
          Mar1[1, i + 1] = Mar1[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
          Mar1[2, i + 1] = Mar1[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
        }
        
      }
      
    }
    
    
  }
  return(list(Mar1 = Mar1, Mar2 = Mar2, Mar3 = Mar3, P = P,Q = Q,  R = R, S = S, U = U, V = V))
}


MARpop = function(P, Q, R, S, alpha, beta, gamma, delta, N, xdebut, ydebut){
  Mar = matrix(0, nrow = 2, ncol = N)
  Mar[1,1] = xdebut
  Mar[2,1] = ydebut
  K = nrow(P)
  for(i in 1:(N-1)){
    U = runif(1)                                                  #On tire un nombre aleatoire entre 0 et 1
    a = P[Mar[1,i], Mar[2,i]]                                       #On créer des variables pour simplifier l'écriture
    b = P[Mar[1,i], Mar[2,i] + K]
    c = Q[Mar[1,i], Mar[2,i] + K]
    d = Q[Mar[1,i], Mar[2,i]]
    if((a > 10**300) | (b > 10**300) | (c > 10**300) | (d > 10**300)){
      warning('renforcement trop important')
      return(list(Mar = matrix(Mar[Mar > 0], nrow = 2, byrow = FALSE), P = P, Q = Q, R = R, S = S, i = i))
    }
    if(U < (a / (a + b + c + d))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
      P[Mar[1,i], Mar[2,i]] = beta * P[Mar[1,i], Mar[2,i]] + alpha                           #On va a gauche et on ajoute delta 
      R[Mar[1,i], Mar[2,i]] = max(0, delta * R[Mar[1,i], Mar[2,i]] - gamma)
      Mar[1, i + 1] = Mar[1, i] - 1
      Mar[2, i + 1] = Mar[2, i]
    }
    else if(U < ((a + b) / (a + b + c + d))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
      P[Mar[1,i], Mar[2,i] + K] = beta * P[Mar[1,i], Mar[2,i] + K] + alpha                                               #On va a droite et on enleve delta
      R[Mar[1,i], Mar[2,i] + K] = max(0, delta * R[Mar[1,i], Mar[2,i] + K] - gamma)
      Mar[1, i + 1] = Mar[1, i] + 1
      Mar[2, i + 1] = Mar[2, i]
    }
    else if(U < ((a + b + c) / (a + b + c + d))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
      Q[Mar[1,i], Mar[2,i] + K] = beta * Q[Mar[1,i], Mar[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
      S[Mar[1,i], Mar[2,i] + K] = max(0, delta * S[Mar[1,i], Mar[2,i] + K] - gamma)
      Mar[1, i + 1] = Mar[1, i]
      Mar[2, i + 1] = Mar[2, i] + 1
    }
    else if((U != 1) | (Mar[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
      Q[Mar[1,i], Mar[2,i]] = beta * Q[Mar[1,i], Mar[2,i]] + alpha                                                      #On va derriere et on enleve delta
      S[Mar[1,i], Mar[2,i]] = max(0, delta * S[Mar[1,i], Mar[2,i]] - gamma)
      Mar[1, i + 1] = Mar[1, i]
      Mar[2, i + 1] = Mar[2, i] - 1
    }
    else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
      Mar[1, i + 1] = Mar[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
      Mar[2, i + 1] = Mar[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
    }
  }
  return(list(Mar = Mar, P = P, Q = Q,  R = R, S = S))
}


MARpopBis = function(P, Q, R, S, U, V, alpha, beta, gamma, delta, N, xdebut, ydebut){
  Mar = matrix(0, nrow = 2, ncol = N)
  Mar[1,1] = xdebut
  Mar[2,1] = ydebut
  K = nrow(P)
  for(i in 1:(N-1)){
    u = runif(1)                                                  #On tire un nombre aleatoire entre 0 et 1
    a = P[Mar[1,i], Mar[2,i]]                                       #On créer des variables pour simplifier l'écriture
    b = P[Mar[1,i], Mar[2,i] + K]
    c = Q[Mar[1,i], Mar[2,i] + K]
    d = Q[Mar[1,i], Mar[2,i]]
    if((a > 10**300) | (b > 10**300) | (c > 10**300) | (d > 10**300)){
      warning('renforcement trop important')
      return(list(Mar = matrix(Mar[Mar > 0], nrow = 2, byrow = FALSE), P = P, Q = Q, R = R, S = S, U = U, V = V, i = i))
    }
    if(u < (a / (a + b + c + d))){        #Si U est plus petit que la proba d'aller a gauche (Poid d'aller a gauche diviser par la somme d'aller a droite, a gauche, devant et derriere)
      P[Mar[1,i], Mar[2,i]] = beta * P[Mar[1,i], Mar[2,i]] + alpha                           #On va a gauche et on ajoute delta 
      R[Mar[1,i], Mar[2,i]] = max(0, delta * R[Mar[1,i], Mar[2,i]] - gamma)
      U[Mar[1,i], Mar[2,i]] = max(0, delta * U[Mar[1,i], Mar[2,i]] - gamma)
      Mar[1, i + 1] = Mar[1, i] - 1
      Mar[2, i + 1] = Mar[2, i]
    }
    else if(u < ((a + b) / (a + b + c + d))){    #Si U est entre les probas d'aller a gauche et d'aller a droite ou a gauche (Poid d'aller a gauche + a droite diviser par la somme d'aller a droite, a gauche, devant et derriere, car si arrive ici c'est que U est plus grand que la proba d'aller a gauche)
      P[Mar[1,i], Mar[2,i] + K] = beta * P[Mar[1,i], Mar[2,i] + K] + alpha                                               #On va a droite et on enleve delta
      R[Mar[1,i], Mar[2,i] + K] = max(0, delta * R[Mar[1,i], Mar[2,i] + K] - gamma)
      U[Mar[1,i], Mar[2,i] + K] = max(0, delta * U[Mar[1,i], Mar[2,i] + K] - gamma)
      Mar[1, i + 1] = Mar[1, i] + 1
      Mar[2, i + 1] = Mar[2, i]
    }
    else if(u < ((a + b + c) / (a + b + c + d))){       #Si U est entre les probas d'aller a gauche ou a droite et d'aller a droite ou a gauche ou devant (Poid d'aller a gauche + a droite + devant diviser par la somme d'aller a droite, a gauche, devant et derriere, car si on arrive ici c'est que U est plus grand que la proba d'aller a gauche ou a droite)
      Q[Mar[1,i], Mar[2,i] + K] = beta * Q[Mar[1,i], Mar[2,i] + K] + alpha                                                      #On va devant et on ajoute delta
      S[Mar[1,i], Mar[2,i] + K] = max(0, delta * S[Mar[1,i], Mar[2,i] + K] - gamma)
      V[Mar[1,i], Mar[2,i] + K] = max(0, delta * V[Mar[1,i], Mar[2,i] + K] - gamma)
      Mar[1, i + 1] = Mar[1, i]
      Mar[2, i + 1] = Mar[2, i] + 1
    }
    else if((u != 1) | (Mar[2,i] != 1)){                                                         #Si U est different de 1 ou que nous ne sommes pas tout en bas
      Q[Mar[1,i], Mar[2,i]] = beta * Q[Mar[1,i], Mar[2,i]] + alpha                                                      #On va derriere et on enleve delta
      S[Mar[1,i], Mar[2,i]] = max(0, delta * S[Mar[1,i], Mar[2,i]] - gamma)
      V[Mar[1,i], Mar[2,i]] = max(0, delta * V[Mar[1,i], Mar[2,i]] - gamma)
      Mar[1, i + 1] = Mar[1, i]
      Mar[2, i + 1] = Mar[2, i] - 1
    }
    else{                                                       #Si U est egal a 1 et que nous sommes tout en bas,                                          
      Mar[1, i + 1] = Mar[1, i]                                 #On choisit de rester sur place sans changer les poids pour ne pas fausser l'aleatoire (Juste pour cette iteration).
      Mar[2, i + 1] = Mar[2, i]                                 #Ce cas n'a quasiment aucune chance d'arriver mais il faut tout de meme l'anticiper
    }
  }
  return(list(Mar = Mar, P = P, Q = Q,  R = R, S = S, U = U, V = V))
}


#### Matrices de transitions ####

K = 400                          #On a une grille de taille K*K


P1 = matrix(1,nrow=K,ncol=2)      #Matrice pour les MAR à une dimension
P1[1,1]= 0
P1[K,2]= 0


P = matrix(1,nrow=K,ncol=2*K)    #Matrice des mouvements horizontaux
for(i in 1:K){
  P[1, i] = 0
  P[K, K+i] = 0
}


Q = matrix(1,nrow=K,ncol=2*K)    #Matrice des mouvements verticaux
for(i in 1:K){
  Q[i,1] = 0
  Q[i, 2*K] = 0
}


R = matrix(1,nrow=K,ncol=2*K)    #Matrice des mouvements horizontaux
for(i in 1:K){
  R[1, i] = 0
  R[K, K+i] = 0
}


S = matrix(1,nrow=K,ncol=2*K)    #Matrice des mouvements verticaux
for(i in 1:K){
  S[i,1] = 0
  S[i, 2*K] = 0
}


U = matrix(1,nrow=K,ncol=2*K)    #Matrice des mouvements horizontaux
for(i in 1:K){
  U[1, i] = 0
  U[K, K+i] = 0
}


V = matrix(1,nrow=K,ncol=2*K)    #Matrice des mouvements verticaux
for(i in 1:K){
  V[i,1] = 0
  V[i, 2*K] = 0
}


#### Initialisations et Tests des fonctions ####


debut = floor(K/2) 
alpha = 0.5
beta = 1
delta = 0.999
gamma = 0

MAR(P1, alpha, beta, 100, debut)
MAR2d(P, Q, alpha, beta, 100, debut)
MARri(P, Q, alpha, beta, 100, debut)
MARstop(P, Q, alpha, beta, 100, debut)
MAR2pop(P, Q, R, S, alpha, beta, gamma, delta, 100, debut, debut, debut, debut)
MARpop(P, Q, R, S, alpha, beta, gamma, delta, 100, debut, debut)
MAR3pop(P, Q, R, S, U, V, alpha, beta, gamma, delta, 100, debut, debut, debut, debut, debut, debut)


#### Première modelisations graphiques ####


N = 50000
alpha = 0.01
beta = 1
delta = 0.999
gamma = 0

X = MAR(P1, alpha, beta, N, debut)$Mar
plot(X ,type='s',main="Marche Aléatoire",
     xlab="Temps",
     ylab="Valeurs")


alpha = 0.2
beta = 1
Y = MAR2d(P, Q, alpha, beta, N, debut)$Mar
Y1 = MARri(P, Q, alpha, beta, N, debut)$Mar
Y2 = MARstop(P, Q, alpha, beta, N, debut)$Mar
Z = MAR2pop(P, Q, R, S, alpha, beta, gamma, delta, N, debut, debut, debut, debut)
Y3 = Z$Mar1
Y4 = Z$Mar2
Z1 = MAR3pop(P, Q, R, S, U, V, alpha, beta, gamma, delta, N, debut, debut, debut, debut, debut, debut)
Y5 = Z1$Mar1
Y6 = Z1$Mar2
Y7 = Z1$Mar3


gradient = c("red","yellow","green", "lightblue","darkblue")
walk = 1:ncol(Y)
N1 = ncol(Y2)

x = Y[1,]
y = Y[2,]
df <- data.frame(walk, x, y)

x1 = Y1[1,]
y1 = Y1[2,]
df1 <- data.frame(walk, x1, y1)

x2 = Y2[1,]
y2 = Y2[2,]
df2 <- data.frame(1:N1, x2, y2)

x3 = Y3[1,]
y3 = Y3[2,]
df3 <- data.frame(walk = 1:N, x = x3, y = y3)

x4 = Y4[1,]
y4 = Y4[2,]
df4 <- data.frame(walk = 1:N, x = x4, y = y4)

df5 = rbind(df3, df4)

x5 = Y5[1,]
y5 = Y5[2,]
df6 <- data.frame(walk = 1:N, x = x5, y = y5)

x6 = Y6[1,]
y6 = Y6[2,]
df7 <- data.frame(walk = 1:N, x = x6, y = y6)

x7 = Y7[1,]
y7 = Y7[2,]
df8 <- data.frame(walk = 1:N, x = x7, y = y7)

df9 = rbind(df6, df7, df8)


plot = ggplot(df, aes(x, y, color = 1:ncol(Y))) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal()
#print(plot)                                   #On remarque qu'il est difficile de voir les zones d'aller retour

#plot = plot + geom_density2d(color = "purple")
#print(plot)                                   #Nous avons donc rajoute ici les bassins d'attractions afin de mieux visualiser les zones d'aggregations

#Pour rajouter les bords
plot = plot + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
plot = plot + labs(title="Evolution d'un individu en 2D",x ="Ouest/Est", y = "Sud/Nord")
print(plot)



plot1 = ggplot(df1, aes(x1, y1, color = 1:N)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal()
#print(plot1)                                  #On remarque qu'il est difficile de voir les zones d'aller retour

#plot1 = plot1 + geom_density2d(color = "purple")
#print(plot1)                                  #Nous avons donc rajoute ici les bassins d'attractions afin de mieux visualiser les zones d'aggregations

# Pour rajouter les bords
plot1 = plot1 + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
print(plot1)



plot2 = ggplot(df2, aes(x2, y2, color = 1:N1)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal()
#print(plot1)                                  #On remarque qu'il est difficile de voir les zones d'aller retour

#plot2 = plot1 + geom_density2d(color = "purple")
#print(plot1)                                  #Nous avons donc rajoute ici les bassins d'attractions afin de mieux visualiser les zones d'aggregations

# Pour rajouter les bords
plot2 = plot2 + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
print(plot2)



plot3 = ggplot(df5, aes(x, y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal()
#print(plot1)                                  #On remarque qu'il est difficile de voir les zones d'aller retour

#plot2 = plot1 + geom_density2d(color = "purple")
#print(plot1)                                  #Nous avons donc rajoute ici les bassins d'attractions afin de mieux visualiser les zones d'aggregations

# Pour rajouter les bords
plot3 = plot3 + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
print(plot3)



plot4 = ggplot(df9, aes(x, y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal()
#print(plot4)                                  #On remarque qu'il est difficile de voir les zones d'aller retour

#plot4 = plot4 + geom_density2d(color = "purple")
#print(plot1)                                  #Nous avons donc rajoute ici les bassins d'attractions afin de mieux visualiser les zones d'aggregations

# Pour rajouter les bords
plot4 = plot4 + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
print(plot4)


#### Difference renforcement, sans renforcement ####


### 1 dimension (renforcement additif de 0.05)


par(mfrow = c(1,2))
plot(MAR(P1, alpha, beta, N, debut)$Mar, type = 's', ylim = c(0,K))
plot(MAR(P1, 0, 1, N, debut)$Mar, type = 's', ylim = c(0,K))


### 2 dimensions (renforcement additif de 0.05)

Xrenf = MAR2d(P, Q, alpha, beta, N, debut)$Mar
Xnonrenf = MAR2d(P, Q, 0, 1, N, debut)$Mar

gradient = c("red","yellow","green", "lightblue","darkblue")
walk = 1:N

xrenf = Xrenf[1,]
yrenf = Xrenf[2,]
dfrenf <- data.frame(walk, xrenf, yrenf)

xnonrenf = Xnonrenf[1,]
ynonrenf = Xnonrenf[2,]
dfnonrenf <- data.frame(walk, xnonrenf, ynonrenf)


plot = ggplot(dfrenf, aes(xrenf, yrenf, color = 1:N)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
plot1 = ggplot(dfnonrenf, aes(xnonrenf, ynonrenf, color = 1:N)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
grid.arrange(plot, plot1, ncol=2)


### 2 dimensions retour interdit (renforcement additif de 0.05)


Yrenf = MARri(P, Q, alpha, beta, N, debut)$Mar
Ynonrenf = MARri(P, Q, 0, 1, N, debut)$Mar

gradient = c("red","yellow","green", "lightblue","darkblue")
walk = 1:N

x1renf = Yrenf[1,]
y1renf = Yrenf[2,]
df1renf <- data.frame(walk, x1renf, y1renf)

x1nonrenf = Ynonrenf[1,]
y1nonrenf = Ynonrenf[2,]
df1nonrenf <- data.frame(walk, x1nonrenf, y1nonrenf)


plot = ggplot(df1renf, aes(x1renf, y1renf, color = 1:N)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
plot1 = ggplot(df1nonrenf, aes(x1nonrenf, y1nonrenf, color = 1:N)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
grid.arrange(plot, plot1, ncol=2)


### 2 dimensions, arret aux bords (renforcement additif de 0.05)


Zrenf = MARstop(P, Q, alpha, beta, N, debut)$Mar
Znonrenf = MARstop(P, Q, 0, 1, N, debut)$Mar

gradient = c("red","yellow","green", "lightblue","darkblue")
N1 = ncol(Zrenf)
N2 = ncol(Znonrenf)
walk1 = 1:N1
walk2 = 1:N2

x2renf = Zrenf[1,]
y2renf = Zrenf[2,]
df2renf <- data.frame(walk1, x2renf, y2renf)

x2nonrenf = Znonrenf[1,]
y2nonrenf = Znonrenf[2,]
df2nonrenf <- data.frame(walk2, x2nonrenf, y2nonrenf)


plot = ggplot(df2renf, aes(x2renf, y2renf, color = 1:N1)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
plot1 = ggplot(df2nonrenf, aes(x2nonrenf, y2nonrenf, color = 1:N2)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
grid.arrange(plot, plot1, ncol=2)


### Comparaisons des sensibilites aux renforcements
alpha = 3
beta = 1

Z1 = MAR2d(P, Q, alpha, beta, N, debut)$Mar
Z2 = MARri(P, Q, alpha, beta, N, debut)$Mar
Z3 = MARstop(P, Q, alpha, beta, N, debut)$Mar

gradient = c("red","yellow","green", "lightblue","darkblue")
N1 = ncol(Z1)
N2 = ncol(Z2)
N3 = ncol(Z3)
walk1 = 1:N1
walk2 = 1:N2
walk3 = 1:N3

x1 = Z1[1,]
y1 = Z1[2,]
df1 <- data.frame(walk1, x1, y1)

x2 = Z2[1,]
y2 = Z2[2,]
df2 <- data.frame(walk2, x2, y2)

x3 = Z3[1,]
y3 = Z3[2,]
df3 <- data.frame(walk3, x3, y3)


plot1 = ggplot(df1, aes(x1, y1, color = 1:N1)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "MAR2D" ,x="Ouest / Est" , y = "Sud / Nord")
plot2 = ggplot(df2, aes(x2, y2, color = 1:N2)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "MARri" ,x="Ouest / Est" , y = "Sud / Nord")
plot3 = ggplot(df3, aes(x3, y3, color = 1:N3)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "MARstop" ,x="Ouest / Est" , y = "Sud / Nord")
grid.arrange(plot1, plot2, plot3, ncol=3)


#### Variation du renforcement ####


#Renforcement additif :


PetitAlpha = 0.01
GrandAlpha = 0.5


par(mfrow = c(1,2))
plot(MAR(P1, PetitAlpha, 1, N, debut)$Mar, type = 's', ylim = c(0,K))
plot(MAR(P1, GrandAlpha, 1, N, debut)$Mar, type = 's', ylim = c(0,K))


T1 = MAR2d(P, Q, PetitAlpha, 1, N, debut)$Mar
T2 = MAR2d(P, Q, GrandAlpha, 1, N, debut)$Mar

gradient = c("red","yellow","green", "lightblue","darkblue")
walk1 = 1:(ncol(T1))
walk2 = 1:(ncol(T2))

u1 = T1[1,]
v1 = T1[2,]
df1 <- data.frame(walk1, u1, v1)

u2 = T2[1,]
v2 = T2[2,]
df2 <- data.frame(walk2, u2, v2)


plot1 = ggplot(df1, aes(u1, v1, color = 1:ncol(T1))) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
plot2 = ggplot(df2, aes(u2, v2, color = 1:ncol(T2))) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
grid.arrange(plot1, plot2, ncol=2)


#Renforcement multiplicatif : 


PetitBeta = 1.0001
GrandBeta = 1.01

T3 = MAR2d(P, Q, 0, PetitBeta, N, debut)$Mar
T4 = MAR2d(P, Q, 0, GrandBeta, N, debut)$Mar

gradient = c("red","yellow","green", "lightblue","darkblue")
walk3 = 1:(ncol(T3))
walk4 = 1:(ncol(T4))

u3 = T3[1,]
v3 = T3[2,]
df3 <- data.frame(walk3, u3, v3)

u4 = T4[1,]
v4 = T4[2,]
df4 <- data.frame(walk4, u4, v4)


plot3 = ggplot(df3, aes(u3, v3, color = 1:ncol(T3))) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
plot4 = ggplot(df4, aes(u4, v4, color = 1:ncol(T4))) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
grid.arrange(plot3, plot4, ncol=2)


#Renforcement "lineaire"

PetitAlpha = 0.01
GrandAlpha = 0.5
PetitBeta = 1.0001
GrandBeta = 1.01

T5 = MAR2d(P, Q, PetitAlpha, PetitBeta, N, debut)$Mar
T6 = MAR2d(P, Q, GrandAlpha, GrandBeta, N, debut)$Mar

gradient = c("red","yellow","green", "lightblue","darkblue")
walk5 = 1:(ncol(T5))
walk6 = 1:(ncol(T6))

u5 = T5[1,]
v5 = T5[2,]
df5 <- data.frame(walk5, u5, v5)

u6 = T6[1,]
v6 = T6[2,]
df6 <- data.frame(walk6, u6, v6)


plot5 = ggplot(df5, aes(u5, v5, color = 1:ncol(T5))) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
plot6 = ggplot(df6, aes(u6, v6, color = 1:ncol(T6))) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
grid.arrange(plot5, plot6, ncol=2)


#### Variation de la taille de la grille ####
#Petite grille
K3 = 150                        #On a une grille de taille K*K

P3 = matrix(1,nrow=K3,ncol=2*K3)    #Matrice des mouvements horizontaux
for(i in 1:K3){
  P3[1, i] = 0
  P3[K3, K3+i] = 0
}


Q3 = matrix(1,nrow=K3,ncol=2*K3)    #Matrice des mouvements verticaux
for(i in 1:K3){
  Q3[i,1] = 0
  Q3[i, 2*K3] = 0
}

#Moyenne -
K4 = 250                        #On a une grille de taille K*K

P4 = matrix(1,nrow=K4,ncol=2*K4)    #Matrice des mouvements horizontaux
for(i in 1:K4){
  P4[1, i] = 0
  P4[K4, K4+i] = 0
}


Q4 = matrix(1,nrow=K4,ncol=2*K4)    #Matrice des mouvements verticaux
for(i in 1:K4){
  Q4[i,1] = 0
  Q4[i, 2*K4] = 0
}

#Moyenne +
K5 = 500                        #On a une grille de taille K*K

P5 = matrix(1,nrow=K5,ncol=2*K5)    #Matrice des mouvements horizontaux
for(i in 1:K5){
  P5[1, i] = 0
  P5[K5, K5+i] = 0
}


Q5 = matrix(1,nrow=K5,ncol=2*K5)    #Matrice des mouvements verticaux
for(i in 1:K5){
  Q5[i,1] = 0
  Q5[i, 2*K5] = 0
}

#Grande grille
K6 = 800                        #On a une grille de taille K*K

P6 = matrix(1,nrow=K6,ncol=2*K6)    #Matrice des mouvements horizontaux
for(i in 1:K6){
  P6[1, i] = 0
  P6[K6, K6+i] = 0
}


Q6 = matrix(1,nrow=K6,ncol=2*K6)    #Matrice des mouvements verticaux
for(i in 1:K6){
  Q6[i,1] = 0
  Q6[i, 2*K6] = 0
}

debut3 = floor(K3/2)
debut4 = floor(K4/2) 
debut5 = floor(K5/2)
debut6 = floor(K6/2)


T7 = MAR2d(P3, Q3, alpha, beta, N, debut3)$Mar
T8 = MAR2d(P4, Q4, alpha, beta, N, debut4)$Mar
T9 = MAR2d(P5, Q5, alpha, beta, N, debut5)$Mar
T10 = MAR2d(P6, Q6, alpha, beta, N, debut6)$Mar

gradient = c("red","yellow","green", "lightblue","darkblue")
walk7 = 1:(ncol(T7))
walk8 = 1:(ncol(T8))
walk9 = 1:(ncol(T9))
walk10 = 1:(ncol(T10))

u7 = T7[1,]
v7 = T7[2,]
df7 <- data.frame(walk7, u7, v7)

u8 = T8[1,]
v8 = T8[2,]
df8 <- data.frame(walk8, u8, v8)

u9 = T9[1,]
v9 = T9[2,]
df9 <- data.frame(walk9, u9, v9)

u10 = T10[1,]
v10 = T10[2,]
df10 <- data.frame(walk10, u10, v10)


plot7 = ggplot(df7, aes(u7, v7, color = 1:ncol(T7))) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K3, yend = 0)) + geom_segment(aes(x = 0, y = K3, xend = K3, yend = K3)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K3)) + geom_segment(aes(x = K3, y = 0, xend = K3, yend = K3))+ labs (title = "Petite grille" ,x="Ouest / Est" , y = "Sud / Nord")
plot8 = ggplot(df8, aes(u8, v8, color = 1:ncol(T8))) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K4, yend = 0)) + geom_segment(aes(x = 0, y = K4, xend = K4, yend = K4)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K4)) + geom_segment(aes(x = K4, y = 0, xend = K4, yend = K4))+ labs (title = "Moyenne -" ,x="Ouest / Est" , y = "Sud / Nord")
plot9 = ggplot(df9, aes(u9, v9, color = 1:ncol(T9))) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K5, yend = 0)) + geom_segment(aes(x = 0, y = K5, xend = K5, yend = K5)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K5)) + geom_segment(aes(x = K5, y = 0, xend = K5, yend = K5))+ labs (title = "Moyenne +" ,x="Ouest / Est" , y = "Sud / Nord")
plot10 = ggplot(df10, aes(u10, v10, color = 1:ncol(T10))) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K6, yend = 0)) + geom_segment(aes(x = 0, y = K6, xend = K6, yend = K6)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K6)) + geom_segment(aes(x = K6, y = 0, xend = K6, yend = K6))+ labs (title = "Grande grille" ,x="Ouest / Est" , y = "Sud / Nord")
grid.arrange(plot7, plot8, plot9, plot10, ncol=2)



#### Initialisation pour differents tests ####
# Initialisation
K = 400                          #On a une grille de taille K*K
K1 = 100
debut = floor(K/2)
debut1 = floor(K1/2)

N = 50000
t = seq(2,10000, by = 1000)
t1 = 10000

n = 50
d = 10

alpha = 1
beta = 1
alpha1 = 0
beta1 = 1.01

Alpha = seq(0, 5, by = 0.5)
Beta = seq(1, 1.2, by = 0.01)

# Matrices #
P = matrix(1,nrow=K,ncol=2*K)    #Matrice des mouvements horizontaux
for(i in 1:K){
  P[1, i] = 0
  P[K, K+i] = 0
}
Q = matrix(1,nrow=K,ncol=2*K)    #Matrice des mouvements verticaux
for(i in 1:K){
  Q[i,1] = 0
  Q[i, 2*K] = 0
}
P2 = matrix(1,nrow=K1,ncol=2*K1)    #Matrice des mouvements horizontaux
for(i in 1:K1){
  P2[1, i] = 0
  P2[K1, K1+i] = 0
}
Q2 = matrix(1,nrow=K1,ncol=2*K1)    #Matrice des mouvements verticaux
for(i in 1:K1){
  Q2[i,1] = 0
  Q2[i, 2*K1] = 0
}



ArreteP = array(0, dim = c(K1, 2*K1, n))
ArreteQ = array(0, dim = c(K1, 2*K1, n))

x = 5   #Les x sommets/arretes les plus visités



#### Temps pour toucher le bord ####

#Maximum de N = 50 000 iterations 

#Renforcement additif de 0 a 10 par pas de 0.5 (Tb1), multiplicatif de 1 a 3 par pas de 0.1 (Tb2)


N = 10000
Alpha = seq (0,10, by =0.5 )

Tb1 = rep(NA, length(Alpha))
for(i in 1:length(Alpha)){
  Tb1[i] = mean(replicate(n,dim(MARstop(P2, Q2, Alpha[i], beta, N, debut1)$Mar))[2,])
}

Tb2 = rep(NA, length(Beta))
for(i in 1:length(Beta)){
  Tb2[i] = mean(replicate(n,dim(MARstop(P2, Q2, alpha1, Beta[i], N, debut1)$Mar))[2,])
} 

par(mfrow=c(1,1))
plot(Alpha, Tb1, type = 'l',ylab="Temps pour toucher le bord",main="Temps pour toucher le bord en fonction du alpha")
plot(Beta, Tb2, type = 'l',ylab="Temps pour toucher le bord",main="Temps pour toucher le bord en fonction du beta")       #On remarque que certaines iterations n'ont pas touche le bord mais se sont ronforce trop vite. Ce graphe n'a donc pas vraiment de sens car il montre en combien d'iteration le renforcement devient trop grand et non le temps pour toucher un bord (qui converge vers l'infini)

#Valeurs de alpha critiques
N = 50000
Alpha = seq (0,15, by = 0.5 )

Tb1 = rep(NA, length(Alpha))
for(i in 1:length(Alpha)){
  Tb1[i] = mean(replicate(n,dim(MARstop(P2, Q2, Alpha[i], beta, 10000, debut1)$Mar))[2,])
}

Tb2 = rep(NA, length(Alpha))
for(i in 1:length(Alpha)){
  Tb2[i] = mean(replicate(n,dim(MARstop(P2, Q2, Alpha[i], beta, 10000, debut1)$Mar))[2,])
} 

par(mfrow=c(1,2))
plot(Alpha, Tb1, type = 'l',ylab="Temps pour toucher le bord",main="Pour N=10.000")
plot(Alpha, Tb2, type = 'l',ylab="Temps pour toucher le bord",main="Pour N=50.000")


#Valeurs de alpha critiques -- Tailles grilles
N = 50000
n=10
Alpha = seq (0,10, by = 0.5 )

Tb1 = rep(NA, length(Alpha))
for(i in 1:length(Alpha)){
  Tb1[i] = mean(replicate(n,dim(MARstop(P2, Q2, Alpha[i], beta, 10000, debut1)$Mar))[2,])
}

Tb2 = rep(NA, length(Alpha))
for(i in 1:length(Alpha)){
  Tb2[i] = mean(replicate(n,dim(MARstop(P, Q, Alpha[i], beta, 10000, debut)$Mar))[2,])
} 

par(mfrow=c(1,2))
plot(Alpha, Tb1, type = 'l',ylab="Temps pour toucher le bord",main="Pour petite grille K=100")
plot(Alpha, Tb2, type = 'l',ylab="Temps pour toucher le bord",main="Pour grande grille K=400")


#### Temps d'arret ####

par(mfrow = c(3,2))
n=50
d=50

#On cherche le nombre d'iteration necessaire pour depasser une distance de d
#MAR2d additive

Ta1 = rep(0, length(Alpha))
for(k in 1:length(Alpha)){
  Tatemp = rep(0,n)  
  X = replicate(n, MAR2d(P, Q, Alpha[k], beta, t1, debut)$Mar)
  for(j in 1:n){  
    i = 1
    while((abs(X[1,i,j] - debut) + abs(X[2,i,j] - debut) < d) & (i < t1)){
      i = i + 1
    }
    if(i >= t1){
      Tatemp[j] = 0
    }
    else{
      Tatemp[j] = i
    }
  }
  Ta1[k] = mean(Tatemp)
}
plot(Alpha, Ta1, type = 'l',ylab="Temps d'arrêt",main="MAR2d additive et d = 50")


#MARri additive

Ta2 = rep(0, length(Alpha))
for(k in 1:length(Alpha)){
  Tatemp = rep(0,n)  
  X = replicate(n, MARri(P, Q, Alpha[k], beta, t1, debut)$Mar)
  for(j in 1:n){  
    i = 1
    while((abs(X[1,i,j] - debut) + abs(X[2,i,j] - debut) < d) & (i < t1)){
      i = i + 1
    }
    if(i >= t1){
      Tatemp[j] = 0
    }
    else{
      Tatemp[j] = i
    }
  }
  Ta2[k] = mean(Tatemp)
}
plot(Alpha, Ta2, type = 'l',ylab="Temps d'arrêt",main="MARri additive et d = 50")


#### Distance du point de depart apres N iterations ####

par(mfrow = c(2,3))

#MAR2d additive, variation des iterations

dist1 = rep(NA, length(t))
for(i in 1:length(t)){
  X = replicate(n, MAR2d(P, Q, alpha, beta, t[i], debut)$Mar)
  distance = sqrt((abs((X[1, t[i],] - debut)))**2 + (abs((X[2, t[i],] - debut)))**2)
  dist1[i] = mean(distance)
}
plot(t, dist1, type = 'l',xlab="Nombre d'itérations",ylab="Distance",main="MAR2D")


#MAR2d additive, variation du renforcement

dist2 = rep(NA, length(Alpha))
for(i in 1:length(Alpha)){
  X = replicate(n, MAR2d(P, Q, Alpha[i], beta, t1, debut)$Mar)
  distance = sqrt((abs((X[1, t1,] - debut)))**2 + (abs((X[2, t1,] - debut)))**2)
  dist2[i] = mean(distance)
}
plot(Alpha, dist2, type = 'l',xlab="Alpha",ylab="Distance",main="MAR2D")

#MAR2d non renforcé, variation des iterations

dist3 = rep(NA, length(t))
for(i in 1:length(t)){
  X = replicate(n, MAR2d(P, Q, alpha1, beta, t[i], debut)$Mar)
  distance = sqrt((abs((X[1, t[i],] - debut)))**2 + (abs((X[2, t[i],] - debut)))**2)
  dist3[i] = mean(distance)
}
plot(t, dist3, type = 'l',xlab="Nombre d'itérations",ylab="Distance",main="Marche Aléatoire 2d")

#MARri additive, variation des iterations

dist4 = rep(NA, length(t))
for(i in 1:length(t)){
  X = replicate(n, MARri(P, Q, alpha, beta, t[i], debut)$Mar)
  distance = sqrt((abs((X[1, t[i],] - debut)))**2 + (abs((X[2, t[i],] - debut)))**2)
  dist4[i] = mean(distance)
}
plot(t, dist4, type = 'l',xlab="Nombre d'itérations",ylab="Distance",main="MARri")

#MARri additive, variation du renforcement

dist5 = rep(NA, length(Alpha))
for(i in 1:length(Alpha)){
  X = replicate(n, MARri(P, Q, Alpha[i], beta, t1, debut)$Mar)
  distance = sqrt((abs((X[1, t1,] - debut)))**2 + (abs((X[2, t1,] - debut)))**2)
  dist5[i] = mean(distance)
}
plot(Alpha, dist5, type = 'l',xlab="Alpha",ylab="Distance",main="MARri")

#MARri non renforcé, variation des iterations

dist6 = rep(NA, length(t))
for(i in 1:length(t)){
  X = replicate(n, MAR2d(P, Q, alpha1, beta, t[i], debut)$Mar)
  distance = sqrt((abs((X[1, t[i],] - debut)))**2 + (abs((X[2, t[i],] - debut)))**2)
  dist6[i] = mean(distance)
}
plot(t, dist6, type = 'l',xlab="Nombre d'itérations",ylab="Distance",main="Marche Aléatoire RI")


#### Distance max ####

par(mfrow = c(2,3))

#MAR2d additive, variation des iterations

MeanDistMax = rep(0, length(t))
for(l in 1:length(t)){
  X6 = replicate(n, MAR2d(P, Q, alpha, beta, t[l], debut)$Mar)
  DistMax = rep(0,n)
  for(i in 1:n){
    for(j in 1:t[l]){
      Dist = sqrt((abs((X6[1, j, i] - debut)))**2 + (abs((X6[2, j, i] - debut)))**2)
      if(Dist > DistMax[i]){
        DistMax[i] = Dist
        }
      }
    }
  MeanDistMax[l] = mean(DistMax)
  }

plot(t, MeanDistMax, type = 'l')

#MAR2d additive, variation du renforcement

MeanDistMax = rep(0, length(Alpha))
for(l in 1:length(Alpha)){
  X6 = replicate(n, MAR2d(P, Q, Alpha[l], beta1, t1, debut)$Mar)
  DistMax = rep(0,n)
  for(i in 1:n){
    for(j in 1:t1){
      Dist = sqrt((abs((X6[1, j, i] - debut)))**2 + (abs((X6[2, j, i] - debut)))**2)
      if(Dist > DistMax[i]){
        DistMax[i] = Dist
      }
    }
  }
  MeanDistMax[l] = mean(DistMax)
}

plot(Alpha, MeanDistMax, type = 'l')

#MAR2d non renforcé, variation des iterations

MeanDistMax = rep(0, length(t))
for(l in 1:length(t)){
  X6 = replicate(n, MAR2d(P, Q, alpha1, beta, t[l], debut)$Mar)
  DistMax = rep(0,n)
  for(i in 1:n){
    for(j in 1:t[l]){
      Dist = sqrt((abs((X6[1, j, i] - debut)))**2 + (abs((X6[2, j, i] - debut)))**2)
      if(Dist > DistMax[i]){
        DistMax[i] = Dist
      }
    }
  }
  MeanDistMax[l] = mean(DistMax)
}

plot(t, MeanDistMax, type = 's')

#MARri additive, variation des iterations

MeanDistMax = rep(0, length(t))
for(l in 1:length(t)){
  X6 = replicate(n, MARri(P, Q, alpha, beta, t[l], debut)$Mar)
  DistMax = rep(0,n)
  for(i in 1:n){
    for(j in 1:t[l]){
      Dist = sqrt((abs((X6[1, j, i] - debut)))**2 + (abs((X6[2, j, i] - debut)))**2)
      if(Dist > DistMax[i]){
        DistMax[i] = Dist
      }
    }
  }
  MeanDistMax[l] = mean(DistMax)
}

plot(t, MeanDistMax, type = 's')

#MARri additive, variation du renforcement

MeanDistMax = rep(0, length(Alpha))
for(l in 1:length(Alpha)){
  X6 = replicate(n, MARri(P, Q, Alpha[l], beta1, t1, debut)$Mar)
  DistMax = rep(0,n)
  for(i in 1:n){
    for(j in 1:t1){
      Dist = sqrt((abs((X6[1, j, i] - debut)))**2 + (abs((X6[2, j, i] - debut)))**2)
      if(Dist > DistMax[i]){
        DistMax[i] = Dist
      }
    }
  }
  MeanDistMax[l] = mean(DistMax)
}

plot(Alpha, MeanDistMax, type = 'l')

#MARri non renforcé, variation des iterations

MeanDistMax = rep(0, length(t))
for(l in 1:length(t)){
  X6 = replicate(n, MARri(P, Q, alpha1, beta, t[l], debut)$Mar)
  DistMax = rep(0,n)
  for(i in 1:n){
    for(j in 1:t[l]){
      Dist = sqrt((abs((X6[1, j, i] - debut)))**2 + (abs((X6[2, j, i] - debut)))**2)
      if(Dist > DistMax[i]){
        DistMax[i] = Dist
      }
    }
  }
  MeanDistMax[l] = mean(DistMax)
}

plot(t, MeanDistMax, type = 's')

#### Aretes les plus empruntees ####


#On choisit quel type de MAR on fait ici (pour les arretes les plus empruntées) :


for(i in 1:n){
  A = MAR2d(P2, Q2, alpha, beta, 100, debut1)
  ArreteP[1:K1, 1:(2*K1),i] = A$ArreteP
  ArreteQ[1:K1, 1:(2*K1),i] = A$ArreteQ
}

#On commence par regarder les indices des aretes les plus empruntees pour chaque simulation

idP = matrix(0, nrow = n, ncol = 2)
idQ = matrix(0, nrow = n, ncol = 2)


for(i in 1:n){
  idp = which.max(ArreteP[1:K1, 1:(2*K1), i])
  if((idp / K1) - (floor(idp / K1)) == 0){
    idP[i, 1] = K1
    idP[i, 2] = idp / K1
  } else{
    idP[i, 1] = ((idp / K1) - floor(idp / K1))*K1
    idP[i, 2] = floor(idp / K1) + 1
  }
  idq = which.max(ArreteQ[1:K1, 1:(2*K1), i])
  if((idq / K1) - (floor(idq / K1)) == 0){
    idQ[i, 1] = K1
    idQ[i, 2] = idq / K1
  } else{
    idQ[i, 1] = ((idq / K1) - floor(idq / K1))*K1
    idQ[i, 2] = floor(idq / K1) + 1
  }
}

idP     #On obtient les sommets par lesquels partent chaque arrête. Si la deuxième colone est > K1, on va a droite, sinon on va a gauche
idQ     #On obtient les sommets par lesquels partent chaque arrête. Si la deuxième colone est > K1, on monte, sinon on descend

#On regarde maintenant les aretes les plus empruntees toutes simulation confondu. Pour cela on commence par sommer les matrices de chaque simulations

arreteP = matrix(0, nrow = K1, ncol = 2*K1)
arreteQ = matrix(0, nrow = K1, ncol = 2*K1)

for(i in 1:K1){
  for(j in 1:K1){
    arreteP[i, j] = sum(ArreteP[i, j, 1:n]) 
    arreteP[i, j + K1] = sum(ArreteP[i, j + K1, 1:n])
    arreteQ[i, j] = sum(ArreteQ[i, j, 1:n])
    arreteQ[i, j + K1] = sum(ArreteQ[i, j + K1, 1:n])
  }
}

#On cherche maintenant les x aretes les plus empruntees (horizontalement et verticalement)
#On va d'abord copie les matrices crees precedemment et a chaque fois qu'on obtient l'indice du max, on remplace la valeur du max par 0 dans la copie
#Pour obtenir au final les x aretes les plus empruntees

CarreteP = arreteP
CarreteQ = arreteQ

id1P = matrix(0, nrow = x, ncol = 4)
id1Q = matrix(0, nrow = x, ncol = 4)

for(i in 1:x){
  id1p = which.max(CarreteP)
  id1P[i, 1] = (id1p - 1)%%nrow(CarreteP) + 1
  id1P[i, 2] = ceiling(id1p/nrow(CarreteP))
  CarreteP[id1P[i, 1], id1P[i, 2]] = 0
  
  id1q = which.max(CarreteQ)
  id1Q[i, 1] = (id1q - 1)%%nrow(CarreteQ) + 1
  id1Q[i, 2] = ceiling(id1q/nrow(CarreteQ))
  CarreteQ[id1Q[i, 1], id1Q[i, 2]] = 0
  
  if(id1P[i, 2] > K1){
    id1P[i, 2] = id1P[i, 2] - K1
    id1P[i, 3] = 2
    id1P[i, 4] = arreteP[id1P[i,1], id1P[i,2] + K1]
  }
  else{
    id1P[i, 3] = 1
    id1P[i, 4] = arreteP[id1P[i,1], id1P[i,2]]
  }
  
  if(id1Q[i, 2] > K1){
    id1Q[i, 2] = id1Q[i, 2] - K1
    id1Q[i, 3] = 2
    id1Q[i, 4] = arreteQ[id1Q[i,1], id1Q[i,2] + K1]
  }
  else{
    id1Q[i, 3] = 1
    id1Q[i, 4] = arreteQ[id1Q[i,1], id1Q[i,2]]
  }
}

#On remarque que les aretes les plus empruntees sont proches les unes des autres
id1P    
id1Q


#Plot A
gradient = c("red","yellow","green", "lightblue","darkblue")
walk = 1:ncol(A$Mar)
N = ncol(A$Mar)

x = A[1,]
y = A[2,]
df <- data.frame(walk, x, y)

plot = ggplot(df, aes(x, y, color = 1:ncol(A$Mar))) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal()+ geom_segment(aes(x = 0, y = 0, xend = K1, yend = 0)) + geom_segment(aes(x = 0, y = K1, xend = K1, yend = K1)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K1)) + geom_segment(aes(x = K1, y = 0, xend = K1, yend = K1))+ labs(title="Evolution d'un individu en 2D",x ="Ouest/Est", y = "Sud/Nord")
print(plot)





#### Sommets les plus visités ####


A = replicate(3 ,MAR2d(P, Q, alpha, beta, 200, debut)$Sommet)

Sommet = matrix(0, nrow = x, ncol = 2)
Temp = matrix(0, nrow = K, ncol = K)

for(j in 1:K){
  for(k in 1:K){
    Temp[j,k] = sum(A[j, k, 1:3])
  }
}

for(l in 1:x){
  idMax = which.max(Temp)
  Sommet[l, 1] = (idMax - 1)%%nrow(Temp) + 1
  Sommet[l, 2] = ceiling(idMax/nrow(Temp))
  Temp[Sommet[l, 1], Sommet[l, 2]] = 0
}

Sommet

#### Nombres de sommet visité ####
par(mfrow=c(3,2))

#Mar2d additive, variation des iterations

MeanNbS = rep(0, length(t))

for(i in 1:length(t)){
  A = replicate(n,MAR2d(P, Q, alpha, beta, t[i], debut)$Sommet)
  MeanNbS[i] = length(A[A > 0]) / n
}

PS = MeanNbS / (K*K)

plot(t, MeanNbS, type ='l',xlab="Nombre d'itérations",ylab="Nombre moyen de points visités",main="MAR2D")
plot(t, PS, type ='l',xlab="Nombre d'itérations",ylab="Ratio",main="MAR2D")


#Mar2d additive, variation du Alpha

MeanNbS = rep(0, length(Alpha))

for(i in 1:length(Alpha)){
  A = replicate(n,MAR2d(P, Q, Alpha[i], beta, t1, debut)$Sommet)
  MeanNbS[i] = length(A[A > 0]) / n
}
PS = MeanNbS / (K*K)
plot(Alpha, MeanNbS, type ='l',xlab="Alpha",ylab="Nombre moyen de points visités",main="MAR2D")
plot(Alpha, PS, type ='l',xlab="Alpha",ylab="Ratio",main="MAR2D")


#Mar2d non renforcée, variation des itérations


MeanNbS = rep(0, length(t))

for(i in 1:length(t)){
  A = replicate(n,MAR2d(P, Q, alpha1, beta, t[i], debut)$Sommet)
  MeanNbS[i] = length(A[A > 0]) / n
}

PS = MeanNbS / (K*K)

plot(t, MeanNbS, type ='l',xlab="Nombre d'itérations",ylab="Nombre moyen de points visités",main="MA 2D non renforcée")
plot(t, PS, type ='l',xlab="Nombre d'itérations",ylab="Ratio",main="MA 2D non renforcée")


#Marri additive, variation des iterations
par(mfrow=c(3,2))

MeanNbS = rep(0, length(t))

for(i in 1:length(t)){
  A = replicate(n,MARri(P, Q, alpha, beta, t[i], debut)$Sommet)
  MeanNbS[i] = length(A[A > 0]) / n
}

PS = MeanNbS / (K*K)

plot(t, MeanNbS, type ='l',xlab="Nombre d'itérations",ylab="Nombre moyen de points visités",main="MARri")
plot(t, PS, type ='l',xlab="Nombre d'itérations",ylab="Ratio",main="MARri")


#Marri additive, variation du Alpha


MeanNbS = rep(0, length(Alpha))

for(i in 1:length(Alpha)){
  A = replicate(n,MARri(P, Q, Alpha[i], beta, t1, debut)$Sommet)
  MeanNbS[i] = length(A[A > 0]) / n
}

PS = MeanNbS / (K*K)

plot(Alpha, MeanNbS, type ='l',xlab="Alpha",ylab="Nombre moyen de points visités",main="MARri")
plot(Alpha, PS, type ='l',xlab="Alpha",ylab="Ratio",main="MARri")


#Mar2d non renforcée, variation des itérations


MeanNbS = rep(0, length(t))

for(i in 1:length(t)){
  A = replicate(n,MARri(P, Q, alpha1, beta, t[i], debut)$Sommet)
  MeanNbS[i] = length(A[A > 0]) / n
}

PS = MeanNbS / (K*K)

plot(t, MeanNbS, type ='l',xlab="Nombre d'itérations",ylab="Nombre moyen de points visités",main="MARri non renforcée")
plot(t, PS, type ='l',xlab="Nombre d'itérations",ylab="Ratio",main="MARri non renforcée")


#### Initialisation des populations ####


individu = 10
N = 10000
alpha = 1
beta = 1
delta = 0.8
gamma = 0
K = 250                          #On a une grille de taille K*K
debut = floor(K/2)
xdebut = debut
ydebut = debut
xdebut1 = debut
ydebut1 = debut
xdebut2 = debut 
ydebut2 = debut


P = matrix(1,nrow=K,ncol=2*K)    #Matrice des mouvements horizontaux
for(i in 1:K){
  P[1, i] = 0
  P[K, K+i] = 0
}


Q = matrix(1,nrow=K,ncol=2*K)    #Matrice des mouvements verticaux
for(i in 1:K){
  Q[i,1] = 0
  Q[i, 2*K] = 0
}


R = matrix(1,nrow=K,ncol=2*K)    #Matrice des mouvements horizontaux
for(i in 1:K){
  R[1, i] = 0
  R[K, K+i] = 0
}


S = matrix(1,nrow=K,ncol=2*K)    #Matrice des mouvements verticaux
for(i in 1:K){
  S[i,1] = 0
  S[i, 2*K] = 0
}


U = matrix(1,nrow=K,ncol=2*K)    #Matrice des mouvements verticaux
for(i in 1:K){
  U[1, i] = 0
  U[K, K + i] = 0
}


V = matrix(1,nrow=K,ncol=2*K)    #Matrice des mouvements verticaux
for(i in 1:K){
  V[i,1] = 0
  V[i, 2*K] = 0
}


gradient = c("red","yellow","green", "lightblue","darkblue")
gradient1 = c("orange","red")
gradient2 = c("lightgreen","darkgreen")
gradient3 = c("blue","purple")


X3 = array(0, dim = c(2, N, individu))
X4 = array(0, dim = c(2, N, individu))
X5 = array(0, dim = c(2, N, individu))
X15 = array(0, dim = c(2, N, individu))
X16 = array(0, dim = c(2, N, individu))
X17 = array(0, dim = c(2, N, individu))

X39 = array(0, dim = c(2, N, individu))
X40 = array(0, dim = c(2, N, individu))
X47 = array(0, dim = c(2, N, individu))
X48 = array(0, dim = c(2, N, individu))
X49 = array(0, dim = c(2, N, individu))

X59 = array(0, dim = c(2, N, individu))
X60 = array(0, dim = c(2, N, individu))
X67 = array(0, dim = c(2, N, individu))
X68 = array(0, dim = c(2, N, individu))
X69 = array(0, dim = c(2, N, individu))


#### 1 Population ####


for(i in 1:(individu)){
  A = MAR2d(P, Q, alpha, beta, N, debut)
  P = A$P
  Q = A$Q
  X3[1:2, 1:N, i] = A$Mar
}

walk = 1:N
x3 = X3[1, 1:N, 1]
y3 = X3[2, 1:N, 1]
df3 = data.frame(walk = walk, x = x3, y = y3)

for(i in 2:individu){
  xtemp = X3[1, 1:N, i]
  ytemp = X3[2, 1:N, i]
  dftemp = data.frame(walk = walk, x = xtemp, y = ytemp)
  df3 = rbind(df3, dftemp)
}

plot3 = ggplot(df3, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Tous les individus" ,x="Ouest / Est" , y = "Sud / Nord")
#plot3 = plot3 + geom_density2d(color = "purple")
print(plot3)



x3 = X3[1, 1:N, 1]
y3 = X3[2, 1:N, 1]
df6 <- data.frame(walk = walk, x = x3, y = y3)

x3 = X3[1, 1:N, floor(individu/2)]
y3 = X3[2, 1:N, floor(individu/2)]
df7 = data.frame(walk = walk, x = x3, y = y3)

x3 = X3[1, 1:N, individu]
y3 = X3[2, 1:N, individu]
df8 = data.frame(walk = walk, x = x3, y = y3)


plot6 = ggplot(df6, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Premier individu" ,x="Ouest / Est" , y = "Sud / Nord")
plot7 = ggplot(df7, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Individu milieu" ,x="Ouest / Est" , y = "Sud / Nord")
plot8 = ggplot(df8, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Dernier individu" ,x="Ouest / Est" , y = "Sud / Nord")

grid.arrange(plot6, plot7, plot8, plot3, ncol=2, nrow = 2)
grid.arrange(plot6, plot7, plot8, ncol=3)


#### 2 Populations ####


#Pas par pas


for(i in 1:(individu)){
  B = MAR2pop(P, Q, R, S, alpha, beta, gamma, delta, N, xdebut, ydebut, xdebut1, ydebut1)
  P = B$P
  Q = B$Q
  R = B$R
  S = B$S
  X4[1:2, 1:N, i] = B$Mar1
  X5[1:2, 1:N, i] = B$Mar2
}

walk = 1:N
x4 = X4[1, 1:N, 1]
y4 = X4[2, 1:N, 1]
df4 = data.frame(walk = walk, x = x4, y = y4)

x5 = X5[1, 1:N, 1]
y5 = X5[2, 1:N, 1]
df5 = data.frame(walk = walk, x = x5, y = y5)

for(i in 2:individu){
  xtemp1 = X4[1, 1:N, i]
  ytemp1 = X4[2, 1:N, i]
  dftemp1 = data.frame(walk = walk, x = xtemp1, y = ytemp1)
  xtemp2 = X5[1, 1:N, i]
  ytemp2 = X5[2, 1:N, i]
  dftemp2 = data.frame(walk = walk, x = xtemp2, y = ytemp2)
  df4 = rbind(df4, dftemp1)
  df5 = rbind(df5, dftemp2)
}


df4$pop = 1
df5$pop = 2
df79 = rbind(df4,df5)


plot79 = ggplot(df79, aes(x = x, y = y, color = pop)) + geom_path() + labs(color = "Population") + scale_colour_gradientn(colors = c("red","green"))  + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Deux populations" ,x="Ouest / Est" , y = "Sud / Nord")


plot4 = ggplot(df4, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Première population" ,x="Ouest / Est" , y = "Sud / Nord")
#plot4 = plot4 + geom_density2d(color = "purple")

plot5 = ggplot(df5, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Deuxième population" ,x="Ouest / Est" , y = "Sud / Nord")
#plot5 = plot5 + geom_density2d(color = "pink")

grid.arrange(plot4, plot5, ncol=2)
plot79


x4 = X4[1, 1:N, 1]
y4 = X4[2, 1:N, 1]
df9 <- data.frame(walk = walk, x = x4, y = y4)

x4 = X4[1, 1:N, floor(individu/2)]
y4 = X4[2, 1:N, floor(individu/2)]
df10 <- data.frame(walk = walk, x = x4, y = y4)

x4 = X4[1, 1:N, individu]
y4 = X4[2, 1:N, individu]
df11 <- data.frame(walk = walk, x = x4, y = y4)

x5 = X5[1, 1:N, 1]
y5 = X5[2, 1:N, 1]
df12 <- data.frame(walk = walk, x = x5, y = y5)

x5 = X5[1, 1:N, floor(individu/2)]
y5 = X5[2, 1:N, floor(individu/2)]
df13 <- data.frame(walk = walk, x = x5, y = y5)

x5 = X5[1, 1:N, individu]
y5 = X5[2, 1:N, individu]
df14 <- data.frame(walk = walk, x = x5, y = y5)


plot9 = ggplot(df9, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
plot10 = ggplot(df10, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
plot11 = ggplot(df11, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))

plot12 = ggplot(df12, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
plot13 = ggplot(df13, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))
plot14 = ggplot(df14, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))

grid.arrange(plot9, plot10, plot11, plot4, ncol=2, nrow = 2)
grid.arrange(plot9, plot10, plot11, ncol=3)

grid.arrange(plot12, plot13, plot14, plot5, ncol=2, nrow = 2)
grid.arrange(plot12, plot13, plot14, ncol=3)

grid.arrange(plot9, plot10, plot11, plot12, plot13, plot14, ncol=3, nrow = 2)






#Individu par individu





for(i in 1:(individu)){
  B = MARpop(P, Q, R, S, alpha, beta, gamma, delta, N, xdebut, ydebut)
  P = B$P
  Q = B$Q
  R = B$R
  S = B$S
  X59[1:2, 1:N, i] = B$Mar
  B = MARpop(R, S, P, Q, alpha, beta, gamma, delta, N, xdebut1, ydebut1)
  P = B$R
  Q = B$S
  R = B$P
  S = B$Q
  X60[1:2, 1:N, i] = B$Mar
}


walk = 1:N
x59 = X59[1, 1:N, 1]
y59 = X59[2, 1:N, 1]
df59 = data.frame(walk = walk, x = x59, y = y59)

x60 = X60[1, 1:N, 1]
y60 = X60[2, 1:N, 1]
df60 = data.frame(walk = walk, x = x60, y = y60)

for(i in 2:individu){
  xtemp1 = X59[1, 1:N, i]
  ytemp1 = X59[2, 1:N, i]
  dftemp1 = data.frame(walk = walk, x = xtemp1, y = ytemp1)
  xtemp2 = X60[1, 1:N, i]
  ytemp2 = X60[2, 1:N, i]
  dftemp2 = data.frame(walk = walk, x = xtemp2, y = ytemp2)
  df59 = rbind(df59, dftemp1)
  df60 = rbind(df60, dftemp2)
}


df59$pop = 1
df60$pop = 2
df81 = rbind(df59,df60)


plot81 = ggplot(df81, aes(x = x, y = y, color = pop)) + geom_path() + labs(color = "Population") + scale_colour_gradientn(colors = c("red","green"))  + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Deux populations" ,x="Ouest / Est" , y = "Sud / Nord")



plot59 = ggplot(df59, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Première population" ,x="Ouest / Est" , y = "Sud / Nord")
#plot4 = plot4 + geom_density2d(color = "purple")

plot60 = ggplot(df60, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Deuxième population" ,x="Ouest / Est" , y = "Sud / Nord")
#plot5 = plot5 + geom_density2d(color = "purple")

plot81
grid.arrange(plot59, plot60, ncol=2)


x59 = X59[1, 1:N, 1]
y59 = X59[2, 1:N, 1]
df61 <- data.frame(walk = walk, x = x59, y = y59)

x59 = X59[1, 1:N, floor(individu/2)]
y59 = X59[2, 1:N, floor(individu/2)]
df62 <- data.frame(walk = walk, x = x59, y = y59)

x59 = X59[1, 1:N, individu]
y59 = X59[2, 1:N, individu]
df63 <- data.frame(walk = walk, x = x59, y = y59)

x60 = X60[1, 1:N, 1]
y60 = X60[2, 1:N, 1]
df64 <- data.frame(walk = walk, x = x60, y = y60)

x60 = X60[1, 1:N, floor(individu/2)]
y60 = X60[2, 1:N, floor(individu/2)]
df65 <- data.frame(walk = walk, x = x60, y = y60)

x60 = X60[1, 1:N, individu]
y60 = X60[2, 1:N, individu]
df66 <- data.frame(walk = walk, x = x60, y = y60)


plot61 = ggplot(df61, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Premier individu" ,x="Ouest / Est" , y = "Sud / Nord")
plot62 = ggplot(df62, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Milieu" ,x="Ouest / Est" , y = "Sud / Nord")
plot63 = ggplot(df63, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Dernier individu" ,x="Ouest / Est" , y = "Sud / Nord")

plot64 = ggplot(df64, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Premier individu" ,x="Ouest / Est" , y = "Sud / Nord")
plot65 = ggplot(df65, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Milieu" ,x="Ouest / Est" , y = "Sud / Nord")
plot66 = ggplot(df66, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Dernier individu" ,x="Ouest / Est" , y = "Sud / Nord")

grid.arrange(plot61, plot62, plot63, plot59, ncol=2, nrow = 2)
grid.arrange(plot61, plot62, plot63, ncol=3)

grid.arrange(plot64, plot65, plot66, plot60, ncol=2, nrow = 2)
grid.arrange(plot64, plot65, plot66, ncol=3)

grid.arrange(plot61, plot62, plot63, plot64, plot65, plot66, ncol=3, nrow = 2)






#Population par population





for(i in 1:(individu)){
  B = MARpop(P, Q, R, S, alpha, beta, gamma, delta, N, xdebut, ydebut)
  P = B$P
  Q = B$Q
  R = B$R
  S = B$S
  X39[1:2, 1:N, i] = B$Mar
}

for(i in 1:(individu)){
  B = MARpop(R, S, P, Q, alpha, beta, gamma, delta, N, xdebut1, ydebut1)
  P = B$R
  Q = B$S
  R = B$P
  S = B$Q
  X40[1:2, 1:N, i] = B$Mar
}

walk = 1:N
x39 = X39[1, 1:N, 1]
y39 = X39[2, 1:N, 1]
df39 = data.frame(walk = walk, x = x39, y = y39)

x40 = X40[1, 1:N, 1]
y40 = X40[2, 1:N, 1]
df40 = data.frame(walk = walk, x = x40, y = y40)

for(i in 2:individu){
  xtemp1 = X39[1, 1:N, i]
  ytemp1 = X39[2, 1:N, i]
  dftemp1 = data.frame(walk = walk, x = xtemp1, y = ytemp1)
  xtemp2 = X40[1, 1:N, i]
  ytemp2 = X40[2, 1:N, i]
  dftemp2 = data.frame(walk = walk, x = xtemp2, y = ytemp2)
  df39 = rbind(df39, dftemp1)
  df40 = rbind(df40, dftemp2)
}


df39$pop = 1
df40$pop = 2
df80 = rbind(df39,df40)


plot80 = ggplot(df80, aes(x = x, y = y, color = pop)) + geom_path() + labs(color = "Population") + scale_colour_gradientn(colors = c("red","green"))  + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Les deux populations" ,x="Ouest / Est" , y = "Sud / Nord")



plot39 = ggplot(df39, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Première population" ,x="Ouest / Est" , y = "Sud / Nord")
#plot4 = plot4 + geom_density2d(color = "purple")

plot40 = ggplot(df40, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Deuxième population" ,x="Ouest / Est" , y = "Sud / Nord")
#plot5 = plot5 + geom_density2d(color = "purple")

plot80
grid.arrange(plot39, plot40, ncol=2)


x39 = X39[1, 1:N, 1]
y39 = X39[2, 1:N, 1]
df41 <- data.frame(walk = walk, x = x39, y = y39)

x39 = X39[1, 1:N, floor(individu/2)]
y39 = X39[2, 1:N, floor(individu/2)]
df42 <- data.frame(walk = walk, x = x39, y = y39)

x39 = X39[1, 1:N, individu]
y39 = X39[2, 1:N, individu]
df43 <- data.frame(walk = walk, x = x39, y = y39)

x40 = X40[1, 1:N, 1]
y40 = X40[2, 1:N, 1]
df44 <- data.frame(walk = walk, x = x40, y = y40)

x40 = X40[1, 1:N, floor(individu/2)]
y40 = X40[2, 1:N, floor(individu/2)]
df45 <- data.frame(walk = walk, x = x40, y = y40)

x40 = X40[1, 1:N, individu]
y40 = X40[2, 1:N, individu]
df46 <- data.frame(walk = walk, x = x40, y = y40)

#Premiere pop
plot41 = ggplot(df41, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Premier individu" ,x="Ouest / Est" , y = "Sud / Nord")
plot42 = ggplot(df42, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Milieu" ,x="Ouest / Est" , y = "Sud / Nord")
plot43 = ggplot(df43, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Dernier individu" ,x="Ouest / Est" , y = "Sud / Nord")

#Deuxième pop
plot44 = ggplot(df44, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Premier individu" ,x="Ouest / Est" , y = "Sud / Nord")
plot45 = ggplot(df45, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Milieu" ,x="Ouest / Est" , y = "Sud / Nord")
plot46 = ggplot(df46, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Dernier individu" ,x="Ouest / Est" , y = "Sud / Nord")

grid.arrange(plot41, plot42, plot43, plot39, ncol=2, nrow = 2)
grid.arrange(plot41, plot42, plot43, ncol=3)

grid.arrange(plot44, plot45, plot46, plot40, ncol=2, nrow = 2)
grid.arrange(plot44, plot45, plot46, ncol=3)

grid.arrange(plot41, plot42, plot43, plot44, plot45, plot46, ncol=3, nrow = 2)





#### 3 Populations ####





#Pas par pas





for(i in 1:(individu)){
  C = MAR3pop(P, Q, R, S, U, V, alpha, beta, gamma, delta, N, xdebut, ydebut, xdebut1, ydebut1, xdebut2, ydebut2)
  P = C$P
  Q = C$Q
  R = C$R
  S = C$S
  U = C$U
  V = C$V
  X15[1:2, 1:N, i] = C$Mar1
  X16[1:2, 1:N, i] = C$Mar2
  X17[1:2, 1:N, i] = C$Mar3
}

walk = 1:N
x15 = X15[1, 1:N, 1]
y15 = X15[2, 1:N, 1]
df15 = data.frame(walk = walk, x = x15, y = y15)

x16 = X16[1, 1:N, 1]
y16 = X16[2, 1:N, 1]
df16 = data.frame(walk = walk, x = x16, y = y16)

x17 = X17[1, 1:N, 1]
y17 = X17[2, 1:N, 1]
df17 = data.frame(walk = walk, x = x17, y = y17)

for(i in 2:individu){
  xtemp1 = X15[1, 1:N, i]
  ytemp1 = X15[2, 1:N, i]
  dftemp1 = data.frame(walk = walk, x = xtemp1, y = ytemp1)
  xtemp2 = X16[1, 1:N, i]
  ytemp2 = X16[2, 1:N, i]
  dftemp2 = data.frame(walk = walk, x = xtemp2, y = ytemp2)
  xtemp3 = X17[1, 1:N, i]
  ytemp3 = X17[2, 1:N, i]
  dftemp3 = data.frame(walk = walk, x = xtemp3, y = ytemp3)
  df15 = rbind(df15, dftemp1)
  df16 = rbind(df16, dftemp2)
  df17 = rbind(df17, dftemp3)
}


df15$pop = 1
df16$pop = 2
df17$pop = 3
df82 = rbind(df15,df16,df17)


plot82 = ggplot(df82, aes(x = x, y = y, color = pop)) + geom_path() + labs(color = "Population") + scale_colour_gradientn(colors = c("red","green","blue"))  + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Toutes les populations" ,x="Ouest / Est" , y = "Sud / Nord")


plot15 = ggplot(df15, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Première population" ,x="Ouest / Est" , y = "Sud / Nord")
#plot4 = plot4 + geom_density2d(color = "purple")

plot16 = ggplot(df16, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Deuxième population" ,x="Ouest / Est" , y = "Sud / Nord")
#plot5 = plot5 + geom_density2d(color = "purple")

plot17 = ggplot(df17, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Troisième population" ,x="Ouest / Est" , y = "Sud / Nord")
#plot4 = plot4 + geom_density2d(color = "purple")


grid.arrange(plot15, plot16, plot17, ncol=3)
plot82




x15 = X15[1, 1:N, 1]
y15 = X15[2, 1:N, 1]
df18 <- data.frame(walk = walk, x = x15, y = y15)

x15 = X15[1, 1:N, floor(individu/2)]
y15 = X15[2, 1:N, floor(individu/2)]
df19 <- data.frame(walk = walk, x = x15, y = y15)

x15 = X15[1, 1:N, individu]
y15 = X15[2, 1:N, individu]
df20 <- data.frame(walk = walk, x = x15, y = y15)

x16 = X16[1, 1:N, 1]
y16 = X16[2, 1:N, 1]
df21 <- data.frame(walk = walk, x = x16, y = y16)

x16 = X16[1, 1:N, floor(individu/2)]
y16 = X16[2, 1:N, floor(individu/2)]
df22 <- data.frame(walk = walk, x = x16, y = y16)

x16 = X16[1, 1:N, individu]
y16 = X16[2, 1:N, individu]
df23 <- data.frame(walk = walk, x = x16, y = y16)

x17 = X17[1, 1:N, 1]
y17 = X17[2, 1:N, 1]
df24 <- data.frame(walk = walk, x = x17, y = y17)

x17 = X17[1, 1:N, floor(individu/2)]
y17 = X17[2, 1:N, floor(individu/2)]
df25 <- data.frame(walk = walk, x = x17, y = y17)

x17 = X17[1, 1:N, individu]
y17 = X17[2, 1:N, individu]
df26 <- data.frame(walk = walk, x = x17, y = y17)


plot18 = ggplot(df18, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Premier individu PopA" ,x="Ouest / Est" , y = "Sud / Nord")
plot19 = ggplot(df19, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Deuxième individu PopA" ,x="Ouest / Est" , y = "Sud / Nord")
plot20 = ggplot(df20, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Troisième individu PopA" ,x="Ouest / Est" , y = "Sud / Nord")

plot21 = ggplot(df21, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Premier individu PopB" ,x="Ouest / Est" , y = "Sud / Nord")
plot22 = ggplot(df22, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Deuxième individu PopB" ,x="Ouest / Est" , y = "Sud / Nord")
plot23 = ggplot(df23, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Troisième individu PopB" ,x="Ouest / Est" , y = "Sud / Nord")

plot24 = ggplot(df24, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Premier individu PopC" ,x="Ouest / Est" , y = "Sud / Nord")
plot25 = ggplot(df25, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Deuxième individu PopC" ,x="Ouest / Est" , y = "Sud / Nord")
plot26 = ggplot(df26, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Troisième individu PopC" ,x="Ouest / Est" , y = "Sud / Nord")

grid.arrange(plot18, plot19, plot20, plot15, ncol=2, nrow = 2)
grid.arrange(plot18, plot19, plot20, ncol=3)

grid.arrange(plot21, plot22, plot23, plot16, ncol=2, nrow = 2)
grid.arrange(plot21, plot22, plot23, ncol=3)

grid.arrange(plot24, plot25, plot26, plot17, ncol=2, nrow = 2)
grid.arrange(plot24, plot25, plot26, ncol=3)

grid.arrange(plot18, plot19, plot20, plot21, plot22, plot23, plot24, plot25, plot26, ncol=3, nrow = 3)





#Individu par individu





for(i in 1:(individu)){
  C = MARpopBis(P, Q, R, S, U, V, alpha, beta, gamma, delta, N, xdebut, ydebut)
  P = C$P
  Q = C$Q
  R = C$R
  S = C$S
  U = C$U
  V = C$V
  X67[1:2, 1:N, i] = C$Mar
  
  C = MARpopBis(R, S, U, V, P, Q, alpha, beta, gamma, delta, N, xdebut1, ydebut1)
  R = C$P
  S = C$Q
  U = C$R
  V = C$S
  P = C$U
  Q = C$V
  X68[1:2, 1:N, i] = C$Mar
  
  C = MARpopBis(U, V, P, Q, R, S, alpha, beta, gamma, delta, N, xdebut2, ydebut2)
  U = C$P
  V = C$Q
  P = C$R
  Q = C$S
  R = C$U
  S = C$V
  X69[1:2, 1:N, i] = C$Mar
}


walk = 1:N
x67 = X67[1, 1:N, 1]
y67 = X67[2, 1:N, 1]
df67 = data.frame(walk = walk, x = x67, y = y67)

x68 = X68[1, 1:N, 1]
y68 = X68[2, 1:N, 1]
df68 = data.frame(walk = walk, x = x68, y = y68)

x69 = X69[1, 1:N, 1]
y69 = X69[2, 1:N, 1]
df69 = data.frame(walk = walk, x = x69, y = y69)

for(i in 2:individu){
  xtemp1 = X67[1, 1:N, i]
  ytemp1 = X67[2, 1:N, i]
  dftemp1 = data.frame(walk = walk, x = xtemp1, y = ytemp1)
  xtemp2 = X68[1, 1:N, i]
  ytemp2 = X68[2, 1:N, i]
  dftemp2 = data.frame(walk = walk, x = xtemp2, y = ytemp2)
  xtemp3 = X69[1, 1:N, i]
  ytemp3 = X69[2, 1:N, i]
  dftemp3 = data.frame(walk = walk, x = xtemp3, y = ytemp3)
  df67 = rbind(df67, dftemp1)
  df68 = rbind(df68, dftemp2)
  df69 = rbind(df69, dftemp3)
}


df67$pop = 1
df68$pop = 2
df69$pop = 3
df84= rbind(df67,df68,df69)


plot84 = ggplot(df84, aes(x = x, y = y, color = pop)) + geom_path() + labs(color = "Population") + scale_colour_gradientn(colors = c("red","green","blue"))  + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Toutes les populations" ,x="Ouest / Est" , y = "Sud / Nord")


plot67 = ggplot(df67, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Population A" ,x="Ouest / Est" , y = "Sud / Nord")
#plot4 = plot4 + geom_density2d(color = "purple")

plot68 = ggplot(df68, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Population B" ,x="Ouest / Est" , y = "Sud / Nord")
#plot5 = plot5 + geom_density2d(color = "purple")

plot69 = ggplot(df69, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Population C" ,x="Ouest / Est" , y = "Sud / Nord")
#plot4 = plot4 + geom_density2d(color = "purple")

grid.arrange(plot67, plot68, plot69, ncol=3)
plot84


x67 = X67[1, 1:N, 1]
y67 = X67[2, 1:N, 1]
df70 <- data.frame(walk = walk, x = x67, y = y67)

x67 = X67[1, 1:N, floor(individu/2)]
y67 = X67[2, 1:N, floor(individu/2)]
df71 <- data.frame(walk = walk, x = x67, y = y67)

x67 = X67[1, 1:N, individu]
y67 = X67[2, 1:N, individu]
df72 <- data.frame(walk = walk, x = x67, y = y67)

x68 = X68[1, 1:N, 1]
y68 = X68[2, 1:N, 1]
df73 <- data.frame(walk = walk, x = x68, y = y68)

x68 = X68[1, 1:N, floor(individu/2)]
y68 = X68[2, 1:N, floor(individu/2)]
df74 <- data.frame(walk = walk, x = x68, y = y68)

x68 = X68[1, 1:N, individu]
y68 = X68[2, 1:N, individu]
df75 <- data.frame(walk = walk, x = x68, y = y68)

x69 = X69[1, 1:N, 1]
y69 = X69[2, 1:N, 1]
df76 <- data.frame(walk = walk, x = x69, y = y69)

x69 = X69[1, 1:N, floor(individu/2)]
y69 = X69[2, 1:N, floor(individu/2)]
df77 <- data.frame(walk = walk, x = x69, y = y69)

x69 = X69[1, 1:N, individu]
y69 = X69[2, 1:N, individu]
df78 <- data.frame(walk = walk, x = x69, y = y69)


plot70 = ggplot(df70, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Premier Pop A" ,x="Ouest / Est" , y = "Sud / Nord")
plot71 = ggplot(df71, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Milieu Pop A" ,x="Ouest / Est" , y = "Sud / Nord")
plot72 = ggplot(df72, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Dernier Pop A" ,x="Ouest / Est" , y = "Sud / Nord")

plot73 = ggplot(df73, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Premier Pop B" ,x="Ouest / Est" , y = "Sud / Nord")
plot74 = ggplot(df74, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Milieu Pop B" ,x="Ouest / Est" , y = "Sud / Nord")
plot75 = ggplot(df75, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Dernier Pop B" ,x="Ouest / Est" , y = "Sud / Nord")

plot76 = ggplot(df76, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Premier Pop C" ,x="Ouest / Est" , y = "Sud / Nord")
plot77 = ggplot(df77, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Milieu Pop C" ,x="Ouest / Est" , y = "Sud / Nord")
plot78 = ggplot(df78, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Dernier Pop C" ,x="Ouest / Est" , y = "Sud / Nord")

grid.arrange(plot70, plot71, plot72, plot67, ncol=2, nrow = 2)
grid.arrange(plot70, plot71, plot72, ncol=3)

grid.arrange(plot73, plot74, plot75, plot68, ncol=2, nrow = 2)
grid.arrange(plot73, plot74, plot75, ncol=3)

grid.arrange(plot76, plot77, plot78, plot69, ncol=2, nrow = 2)
grid.arrange(plot76, plot77, plot78, ncol=3)

grid.arrange(plot70, plot71, plot72, plot73, plot74, plot75, plot76, plot77, plot78, ncol=3, nrow = 3)





#Population par population





for(i in 1:(individu)){
  C = MARpop(P, Q, R, S, alpha, beta, gamma, delta, N, xdebut, ydebut)
  P = C$P
  Q = C$Q
  R = C$R
  S = C$S
  X47[1:2, 1:N, i] = C$Mar
}

U = C$R
V = C$S

for(i in 1:(individu)){
  C = MARpop(R, S, U, V, alpha, beta, gamma, delta, N, xdebut1, ydebut1)
  R = C$P
  S = C$Q
  U = C$R
  V = C$S
  X48[1:2, 1:N, i] = C$Mar
}

for(i in 1:(individu)){
  C = MARpop(U, V, R, S, alpha, beta, gamma, delta, N, xdebut2, ydebut2)
  U = C$P
  V = C$Q
  X49[1:2, 1:N, i] = C$Mar
}

walk = 1:N
x47 = X47[1, 1:N, 1]
y47 = X47[2, 1:N, 1]
df47 = data.frame(walk = walk, x = x47, y = y47)

x48 = X48[1, 1:N, 1]
y48 = X48[2, 1:N, 1]
df48 = data.frame(walk = walk, x = x48, y = y48)

x49 = X49[1, 1:N, 1]
y49 = X49[2, 1:N, 1]
df49 = data.frame(walk = walk, x = x49, y = y49)

for(i in 2:individu){
  xtemp1 = X47[1, 1:N, i]
  ytemp1 = X47[2, 1:N, i]
  dftemp1 = data.frame(walk = walk, x = xtemp1, y = ytemp1)
  xtemp2 = X48[1, 1:N, i]
  ytemp2 = X48[2, 1:N, i]
  dftemp2 = data.frame(walk = walk, x = xtemp2, y = ytemp2)
  xtemp3 = X49[1, 1:N, i]
  ytemp3 = X49[2, 1:N, i]
  dftemp3 = data.frame(walk = walk, x = xtemp3, y = ytemp3)
  df47 = rbind(df47, dftemp1)
  df48 = rbind(df48, dftemp2)
  df49 = rbind(df49, dftemp3)
}



df47$pop = 1
df48$pop = 2
df49$pop = 3
df83 = rbind(df47,df48,df49)


plot83 = ggplot(df83, aes(x = x, y = y, color = pop)) + geom_path() + labs(color = "Population") + scale_colour_gradientn(colors = c("red","green","blue"))  + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Toutes les populations" ,x="Ouest / Est" , y = "Sud / Nord")




plot47 = ggplot(df47, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Pop A" ,x="Ouest / Est" , y = "Sud / Nord")
#plot4 = plot4 + geom_density2d(color = "purple")

plot48 = ggplot(df48, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Pop B" ,x="Ouest / Est" , y = "Sud / Nord")
#plot5 = plot5 + geom_density2d(color = "purple")

plot49 = ggplot(df49, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Pop C" ,x="Ouest / Est" , y = "Sud / Nord")
#plot4 = plot4 + geom_density2d(color = "purple")

grid.arrange(plot47, plot48, plot49, ncol=3)
plot83




x47 = X47[1, 1:N, 1]
y47 = X47[2, 1:N, 1]
df50 <- data.frame(walk = walk, x = x47, y = y47)

x47 = X47[1, 1:N, floor(individu/2)]
y47 = X47[2, 1:N, floor(individu/2)]
df51 <- data.frame(walk = walk, x = x47, y = y47)

x47 = X47[1, 1:N, individu]
y47 = X47[2, 1:N, individu]
df52 <- data.frame(walk = walk, x = x47, y = y47)

x48 = X48[1, 1:N, 1]
y48 = X48[2, 1:N, 1]
df53 <- data.frame(walk = walk, x = x48, y = y48)

x48 = X48[1, 1:N, floor(individu/2)]
y48 = X48[2, 1:N, floor(individu/2)]
df54 <- data.frame(walk = walk, x = x48, y = y48)

x48 = X48[1, 1:N, individu]
y48 = X48[2, 1:N, individu]
df55 <- data.frame(walk = walk, x = x48, y = y48)

x49 = X49[1, 1:N, 1]
y49 = X49[2, 1:N, 1]
df56 <- data.frame(walk = walk, x = x49, y = y49)

x49 = X49[1, 1:N, floor(individu/2)]
y49 = X49[2, 1:N, floor(individu/2)]
df57 <- data.frame(walk = walk, x = x49, y = y49)

x49 = X49[1, 1:N, individu]
y49 = X49[2, 1:N, individu]
df58 <- data.frame(walk = walk, x = x49, y = y49)


plot50 = ggplot(df50, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Premier Pop A" ,x="Ouest / Est" , y = "Sud / Nord")
plot51 = ggplot(df51, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Milieu Pop A" ,x="Ouest / Est" , y = "Sud / Nord")
plot52 = ggplot(df52, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Dernier Pop A" ,x="Ouest / Est" , y = "Sud / Nord")

plot53 = ggplot(df53, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Premier Pop B" ,x="Ouest / Est" , y = "Sud / Nord")
plot54 = ggplot(df54, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Milieu Pop B" ,x="Ouest / Est" , y = "Sud / Nord")
plot55 = ggplot(df55, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Dernier Pop B" ,x="Ouest / Est" , y = "Sud / Nord")

plot56 = ggplot(df56, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Premier Pop C" ,x="Ouest / Est" , y = "Sud / Nord")
plot57 = ggplot(df57, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Milieu Pop C" ,x="Ouest / Est" , y = "Sud / Nord")
plot58 = ggplot(df58, aes(x = x, y = y, color = walk)) + geom_path() + labs(color = "Time") + scale_colour_gradientn(colors = gradient) + theme_minimal() + geom_segment(aes(x = 0, y = 0, xend = K, yend = 0)) + geom_segment(aes(x = 0, y = K, xend = K, yend = K)) + geom_segment(aes(x = 0, y = 0, xend = 0, yend = K)) + geom_segment(aes(x = K, y = 0, xend = K, yend = K))+ labs (title = "Dernier Pop C" ,x="Ouest / Est" , y = "Sud / Nord")

grid.arrange(plot50, plot51, plot52, plot47, ncol=2, nrow = 2)
grid.arrange(plot50, plot51, plot52, ncol=3)

grid.arrange(plot53, plot54, plot55, plot48, ncol=2, nrow = 2)
grid.arrange(plot53, plot54, plot55, ncol=3)

grid.arrange(plot56, plot57, plot58, plot49, ncol=2, nrow = 2)
grid.arrange(plot56, plot57, plot58, ncol=3)

grid.arrange(plot50, plot51, plot52, plot53, plot54, plot55, plot56, plot57, plot58, ncol=3, nrow = 3)




