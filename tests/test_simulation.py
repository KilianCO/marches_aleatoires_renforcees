"""Tests de la simulation."""

import pytest

from marches import Parametres, simuler, simuler_1d
from marches.simulation import BAS, DROITE, GAUCHE, HAUT, poids_initiaux


def pas_successifs(marche):
    return [
        (x2 - x1, y2 - y1)
        for x1, y1, x2, y2 in zip(marche.x, marche.y, marche.x[1:], marche.y[1:])
    ]


def test_les_aretes_qui_sortent_de_la_grille_ont_un_poids_nul():
    taille = 5
    poids = poids_initiaux(taille)

    def w(x, y, direction):
        return poids[4 * (x * taille + y) + direction]

    for i in range(taille):
        assert w(0, i, GAUCHE) == 0 and w(taille - 1, i, DROITE) == 0
        assert w(i, 0, BAS) == 0 and w(i, taille - 1, HAUT) == 0

    assert w(2, 2, GAUCHE) == w(2, 2, DROITE) == w(2, 2, HAUT) == w(2, 2, BAS) == 1
    assert sum(poids) == 4 * taille * taille - 4 * taille


def test_la_marche_reste_dans_la_grille_et_avance_d_un_voisin_a_la_fois():
    resultat = simuler(taille=8, pas=5_000, alpha=0.5, graine=0)
    marche = resultat.marches[0]

    assert len(marche.x) == marche.pas_effectues + 1 == 5_001
    assert all(0 <= x < 8 for x in marche.x) and all(0 <= y < 8 for y in marche.y)
    assert set(pas_successifs(marche)) <= {(-1, 0), (1, 0), (0, 1), (0, -1)}
    assert marche.arret is None


def test_sans_renforcement_les_poids_ne_changent_pas():
    resultat = simuler(taille=10, pas=2_000, alpha=0.0, beta=1.0, graine=1)

    assert resultat.poids[0] == poids_initiaux(10)


def test_chaque_pas_renforce_exactement_une_arete():
    resultat = simuler(taille=10, pas=2_000, alpha=0.25, beta=1.0, graine=2)

    gain = sum(resultat.poids[0]) - sum(poids_initiaux(10))

    assert gain == pytest.approx(2_000 * 0.25)


def test_le_renforcement_concentre_la_marche():
    def sommets(alpha):
        return sum(
            simuler(taille=120, pas=4_000, alpha=alpha, graine=g).sommets_visites()
            for g in range(12)
        )

    assert sommets(3.0) < 0.7 * sommets(0.0)


def test_retour_interdit():
    marche = simuler(
        taille=12, pas=4_000, alpha=1.0, retour_interdit=True, graine=3,
    ).marches[0]

    pas = pas_successifs(marche)

    assert all(
        (a[0] + b[0], a[1] + b[1]) != (0, 0)
        for a, b in zip(pas, pas[1:])
    )


def test_sans_interdiction_les_retours_existent():
    pas = pas_successifs(simuler(taille=12, pas=500, alpha=1.0, graine=3).marches[0])

    assert any((a[0] + b[0], a[1] + b[1]) == (0, 0) for a, b in zip(pas, pas[1:]))


def test_arret_au_bord():
    resultat = simuler(taille=15, pas=100_000, alpha=0.0, arret_au_bord=True, graine=4)
    marche = resultat.marches[0]

    assert marche.arret == "bord"
    assert marche.pas_effectues < 100_000
    assert 0 in marche.position or 14 in marche.position

    # Aucun sommet du bord n'est atteint avant le dernier pas.
    assert all(0 < x < 14 and 0 < y < 14 for x, y in zip(marche.x[:-1], marche.y[:-1]))


def test_un_renforcement_multiplicatif_fort_arrete_la_marche():
    marche = simuler(taille=50, pas=50_000, alpha=0.0, beta=1.5, graine=5).marches[0]

    assert marche.arret == "renforcement"
    assert marche.pas_effectues < 50_000


def test_les_distances_sont_mesurees_depuis_le_depart():
    marche = simuler(taille=40, pas=3_000, alpha=0.0, graine=6).marches[0]

    distances = [
        ((x - 20) ** 2 + (y - 20) ** 2) ** 0.5
        for x, y in zip(marche.x, marche.y)
    ]

    assert marche.depart == (20, 20)
    assert marche.distance_max == pytest.approx(max(distances))
    assert marche.distance_finale == pytest.approx(distances[-1])


def test_visites_et_proportion():
    resultat = simuler(taille=20, pas=1_000, alpha=0.0, graine=7)
    marche = resultat.marches[0]

    assert sum(resultat.visites[0]) == 1_001
    assert resultat.sommets_visites() == len(set(zip(marche.x, marche.y)))
    assert resultat.proportion_visitee() == resultat.sommets_visites() / 400
    assert resultat.passages()[20 // 2][20 // 2] >= 1


def test_la_meme_graine_redonne_la_meme_marche():
    def jouer():
        marche = simuler(taille=30, pas=2_000, alpha=0.7, graine=8).marches[0]
        return marche.x, marche.y

    assert jouer() == jouer()
    assert jouer() != (
        simuler(taille=30, pas=2_000, alpha=0.7, graine=9).marches[0].x,
        simuler(taille=30, pas=2_000, alpha=0.7, graine=9).marches[0].y,
    )


def test_sans_trajectoires_les_mesures_restent_disponibles():
    complet = simuler(taille=30, pas=2_000, alpha=0.7, graine=10)
    leger = simuler(taille=30, pas=2_000, alpha=0.7, graine=10, trajectoires=False)

    assert leger.marches[0].x == []
    assert leger.marches[0].distance_max == complet.marches[0].distance_max
    assert leger.marches[0].position == complet.marches[0].position
    assert leger.sommets_visites() == complet.sommets_visites()


# ------------------------------------------------------------ populations

@pytest.mark.parametrize("ordre", ["pas", "individu", "population"])
def test_chaque_population_fait_marcher_tous_ses_individus(ordre):
    resultat = simuler(
        taille=20, pas=300, populations=3, individus=4, delta=0.8,
        ordre=ordre, graine=11,
    )

    assert [(m.population, m.individu) for m in resultat.marches] == [
        (population, individu)
        for population in range(3)
        for individu in range(4)
    ]
    assert all(m.pas_effectues == 300 for m in resultat.marches)
    assert len(resultat.poids) == len(resultat.visites) == 3


def test_les_individus_d_une_population_heritent_des_poids():
    resultat = simuler(taille=20, pas=500, alpha=0.5, individus=3, graine=12)

    gain = sum(resultat.poids[0]) - sum(poids_initiaux(20))

    assert gain == pytest.approx(3 * 500 * 0.5)


def test_un_passage_affaiblit_l_arete_pour_les_autres_populations():
    resultat = simuler(
        taille=20, pas=400, alpha=1.0, populations=2, delta=0.5,
        ordre="population", graine=13,
    )
    initiaux = poids_initiaux(20)

    # La population 1 marche en second : personne n'a affaibli ses arêtes
    # après son passage, mais elle a affaibli celles de la population 0.
    affaiblies = [
        w for w, w0 in zip(resultat.poids[0], initiaux) if w < w0
    ]

    assert affaiblies
    assert all(w >= 0 for population in resultat.poids for w in population)


def test_gamma_ne_rend_jamais_un_poids_negatif():
    resultat = simuler(
        taille=10, pas=3_000, alpha=0.2, populations=2, gamma=0.6, delta=1.0,
        ordre="pas", graine=14,
    )

    assert min(min(population) for population in resultat.poids) == 0.0
    assert all(m.pas_effectues == 3_000 for m in resultat.marches)


def test_points_de_depart():
    resultat = simuler(
        taille=30, pas=50, populations=2, departs=[(5, 6), (20, 21)], graine=15,
    )

    assert resultat.marches[0].depart == (5, 6)
    assert (resultat.marches[1].x[0], resultat.marches[1].y[0]) == (20, 21)


@pytest.mark.parametrize("options", [
    {"taille": 2},
    {"pas": 0},
    {"alpha": -1},
    {"beta": 0},
    {"delta": -0.1},
    {"populations": 0},
    {"ordre": "inconnu"},
    {"populations": 2, "departs": [(1, 1)]},
    {"departs": [(999, 1)]},
])
def test_reglages_invalides(options):
    with pytest.raises(ValueError):
        simuler(**options)


def test_parametres_ou_options_mais_pas_les_deux():
    assert simuler(Parametres(taille=10, pas=10, graine=0)).marches[0].pas_effectues == 10

    with pytest.raises(TypeError):
        simuler(Parametres(taille=10, pas=10), alpha=1.0)


# ------------------------------------------------------------ une dimension

def test_marche_1d_reste_sur_la_ligne():
    sortie = simuler_1d(taille=20, pas=5_000, alpha=0.0, graine=0)
    positions = sortie["positions"]

    assert len(positions) == 5_001 and sortie["arret"] is None
    assert min(positions) == 0 and max(positions) == 19
    assert all(abs(b - a) == 1 for a, b in zip(positions, positions[1:]))


def test_marche_1d_multiplicative_s_enferme():
    sortie = simuler_1d(taille=200, pas=200_000, alpha=0.0, beta=1.05, graine=1)

    assert sortie["arret"] == "renforcement"
    assert len(set(sortie["positions"][-200:])) == 2
