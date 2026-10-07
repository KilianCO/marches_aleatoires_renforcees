# Marches aléatoires renforcées

Projet de fin de Master 1 de mathématiques appliquées (Université Claude-Bernard Lyon 1, 2022), réalisé par Kilian Collet et Clément Contamin, sous la responsabilité de Frédérique Bienvenüe.

Une marche aléatoire renforcée préfère les chemins qu'elle a déjà empruntés, comme une fourmi qui suit et dépose des phéromones. Le projet simule ces marches sur une grille et étudie leur comportement : trajectoires, temps pour atteindre le bord, exploration, puis interactions entre plusieurs colonies.

Les simulations sont rejouables et modifiables en ligne sur [kilianco.github.io](https://kilianco.github.io/projets/marches-aleatoires-renforcees/).

```
Projet-2.pdf            le rapport (137 pages)
code.r                  le script R d'origine, conservé tel quel
marches/simulation.py   la simulation, en Python
marches/experiences.py  les expériences du rapport, rejouables
marches/figures.py      tracés avec matplotlib (facultatif)
tests/                  tests
```

## Le modèle

Chaque arête orientée de la grille porte un poids, égal à 1 au départ. La probabilité d'emprunter une arête est son poids divisé par la somme des poids des arêtes qui partent du même sommet. Après chaque passage :

```
poids <- beta * poids + alpha
```

Avec `alpha = 0` et `beta = 1`, c'est une marche aléatoire simple.

Quand plusieurs colonies partagent la grille, chacune a ses propres poids. Un passage renforce l'arête pour sa colonie et l'affaiblit pour les autres :

```
poids <- max(0, delta * poids - gamma)
```

## Utiliser la version Python

Python 3.10 ou plus récent. La simulation n'a besoin d'aucune bibliothèque.

```python
from marches import simuler, balayer, executer

# Une marche renforcée de 20 000 pas sur une grille 200 x 200
resultat = simuler(taille=200, pas=20_000, alpha=0.5)
print(resultat.proportion_visitee(), resultat.marches[0].distance_max)

# Deux colonies qui s'évitent
simuler(populations=2, individus=6, pas=4_000, alpha=1.0, delta=0.8, ordre="pas")

# Une mesure moyennée, en faisant varier le renforcement
balayer("alpha", [0, 1, 2, 4], "distance_max", repetitions=30, taille=200, pas=4_000)

# Une expérience du rapport, avec ses réglages
executer("territoires_trois", graine=1)
```

| Paramètre | Rôle |
|---|---|
| `taille`, `pas` | côté de la grille, nombre de pas de chaque marcheur |
| `alpha`, `beta` | renforcement additif et multiplicatif |
| `gamma`, `delta` | affaiblissement additif et multiplicatif des autres colonies |
| `retour_interdit` | le marcheur ne revient pas immédiatement sur son dernier pas |
| `arret_au_bord` | la marche s'arrête en touchant un bord (sinon les bords sont réfléchissants) |
| `populations`, `individus` | nombre de colonies, marcheurs successifs par colonie |
| `ordre` | `"pas"`, `"individu"` ou `"population"` : dans quel ordre les marcheurs avancent |
| `departs` | point de départ de chaque colonie (le centre par défaut) |
| `graine` | graine du générateur aléatoire, pour rejouer la même marche |

Tracer une expérience (demande matplotlib) :

```bash
python -m marches.figures                    # liste les expériences
python -m marches.figures grille_additive
python -m marches.figures territoires_trois --sortie territoires.png
```

Lancer les tests :

```bash
pip install pytest
python -m pytest
```

## Du script R au module Python

Le script R compte 3 000 lignes : chaque variante (`MAR`, `MAR2d`, `MARri`, `MARstop`, `MARpop`, `MAR2pop`, `MAR3pop`) est une fonction écrite séparément, et chaque figure a son propre bloc de code.

La version Python garde le même modèle et le réorganise :

| Dans le script R | En Python |
|---|---|
| sept fonctions de simulation | une fonction `simuler`, dont les variantes sont des paramètres |
| une, deux ou trois populations, codées à la main | un nombre quelconque de colonies |
| matrices `P`, `Q`, `R`, `S`, `U`, `V` de taille K × 2K | un tableau de poids par colonie, quatre directions par sommet |
| un bloc de code par étude statistique | une fonction `balayer` : un paramètre à faire varier, une mesure à moyenner |
| figures produites une à une | un catalogue `EXPERIENCES`, une entrée par figure du rapport |
| pas de graine | simulations reproductibles par une graine |

Deux écarts volontaires avec le script R : une marche dont toutes les arêtes permises ont été affaiblies jusqu'à zéro choisit un voisin au hasard (le script R divisait alors par zéro), et le temps pour toucher le bord est accompagné de la part des marches qui l'ont réellement atteint.

## Principales conclusions du rapport

- Le renforcement comprime la trajectoire : plus `alpha` est grand, moins la marche explore et plus elle met de temps à atteindre le bord.
- Le renforcement multiplicatif s'emballe : dès `beta` proche de 1,15, la marche finit souvent enfermée dans un aller-retour.
- Interdire le retour immédiat fait explorer davantage et rend la marche moins sensible au renforcement.
- On observe des zones d'agrégation plutôt que des routes ; avec plusieurs colonies, ces zones deviennent des territoires.
