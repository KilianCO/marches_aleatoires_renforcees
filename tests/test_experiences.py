"""Tests du catalogue d'expériences."""

import json

import pytest

from marches import EXPERIENCES, balayer, catalogue, executer, mesurer
from marches.simulation import Parametres


# Réglages réduits pour que chaque expérience se joue en un instant.
PETIT = {"taille": 30, "pas": 200}


@pytest.mark.parametrize("identifiant", list(EXPERIENCES))
def test_chaque_experience_se_joue_et_se_convertit_en_json(identifiant):
    experience = EXPERIENCES[identifiant]

    surcharges = dict(PETIT)

    if "departs" in experience["parametres"]:
        surcharges["departs"] = [[8, 8], [20, 20], [8, 20]][: experience["parametres"]["populations"]]

    if experience.get("variable") == "pas":
        del surcharges["pas"]

    sortie = executer(identifiant, surcharges, repetitions=2, graine=0)

    assert sortie["nature"] == experience["nature"]
    assert json.loads(json.dumps(sortie)) == sortie

    if sortie["nature"] == "ligne":
        assert len(sortie["positions"]) == 201

    elif sortie["nature"] == "trajectoire":
        attendu = experience["parametres"].get("populations", 1) * experience["parametres"].get("individus", 1)

        assert len(sortie["marches"]) == attendu
        assert 0 < sortie["proportion_visitee"] <= 1

    else:
        assert [s["nom"] for s in sortie["series"]] == [s["nom"] for s in experience["series"]]
        assert all(len(s["moyennes"]) == len(experience["valeurs"]) for s in sortie["series"])


def test_le_catalogue_decrit_chaque_experience():
    entrees = catalogue()

    assert [e["id"] for e in entrees] == list(EXPERIENCES)
    assert all(e["titre"] and e["figure"] and e["description"] for e in entrees)
    assert {e["nature"] for e in entrees} == {"ligne", "trajectoire", "courbe"}
    assert json.dumps(entrees)


def test_experience_inconnue():
    with pytest.raises(ValueError):
        executer("inconnue")


def test_les_surcharges_remplacent_les_reglages():
    sortie = executer("grille_additive", {"taille": 25, "pas": 40, "alpha": 2.5}, graine=1)

    assert sortie["parametres"]["alpha"] == 2.5
    assert sortie["marches"][0]["pas"] == 40
    assert max(sortie["marches"][0]["x"]) < 25


def test_la_meme_graine_redonne_la_meme_experience():
    assert executer("deux_populations_pas", PETIT, graine=3) == executer("deux_populations_pas", PETIT, graine=3)
    assert executer("temps_bord", PETIT, repetitions=3, graine=3) == executer("temps_bord", PETIT, repetitions=3, graine=3)


def test_balayer_moyenne_la_mesure_sur_les_repetitions():
    sortie = balayer(
        "alpha", [0.0, 4.0], "distance_max",
        repetitions=30, graine=0, taille=120, pas=2_000,
    )

    assert sortie["valeurs"] == [0.0, 4.0] and sortie["repetitions"] == 30
    assert len(sortie["ecarts_types"]) == 2

    # Conclusion du rapport : le renforcement comprime la trajectoire.
    assert sortie["moyennes"][1] < sortie["moyennes"][0]


def test_temps_bord_indique_la_part_des_marches_arrivees():
    sortie = balayer(
        "alpha", [0.0, 50.0], "temps_bord",
        repetitions=10, graine=0, taille=21, pas=3_000, arret_au_bord=True,
    )

    assert sortie["atteint"][0] == 1.0
    assert sortie["moyennes"][0] < sortie["moyennes"][1]
    assert sortie["atteint"][1] < 1.0


def test_le_retour_interdit_fait_explorer_davantage():
    def sommets(retour_interdit):
        return balayer(
            "pas", [3_000], "sommets_visites", repetitions=20, graine=0,
            taille=150, alpha=1.0, retour_interdit=retour_interdit,
        )["moyennes"][0]

    assert sommets(True) > sommets(False)


def test_progres_est_signale():
    fractions = []

    executer("distance_max", PETIT, repetitions=1, graine=0, progres=fractions.append)

    assert fractions == sorted(fractions) and fractions[-1] == pytest.approx(1.0)
    assert len(fractions) == 2 * len(EXPERIENCES["distance_max"]["valeurs"])


def test_mesurer():
    p = Parametres(taille=30, pas=300, alpha=0.5, graine=0, trajectoires=False)

    assert mesurer("sommets_visites", p) >= 1
    assert 0 < mesurer("proportion_visitee", p) <= 1

    with pytest.raises(ValueError):
        mesurer("inconnue", p)


def test_normaliser_accepte_toutes_les_sorties():
    from marches import normaliser, simuler, simuler_1d

    trajectoire = normaliser(simuler(taille=20, pas=50, populations=2, graine=0))
    assert trajectoire["nature"] == "trajectoire" and len(trajectoire["marches"]) == 2
    assert trajectoire["parametres"]["populations"] == 2

    assert normaliser(simuler_1d(taille=20, pas=50, graine=0))["nature"] == "ligne"

    courbe = balayer("alpha", [0, 1], "distance_max", repetitions=2, taille=20, pas=50)
    assert normaliser(courbe)["series"][0]["nom"] == "Série 1"
    assert len(normaliser([courbe, {**courbe, "nom": "B"}])["series"]) == 2

    deja = executer("ligne_simple", PETIT, graine=0)
    assert normaliser(deja) is deja

    with pytest.raises(TypeError):
        normaliser(42)


def test_les_figures_se_tracent(tmp_path):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    import matplotlib.pyplot as plt

    from marches.figures import tracer

    for identifiant in ("ligne_additive", "territoires_deux", "distance_max"):
        surcharges = dict(PETIT)

        if identifiant == "territoires_deux":
            surcharges["departs"] = [[8, 8], [20, 20]]

        axe = tracer(executer(identifiant, surcharges, repetitions=2, graine=0))
        axe.figure.savefig(tmp_path / f"{identifiant}.png")
        plt.close(axe.figure)

        assert (tmp_path / f"{identifiant}.png").stat().st_size > 2_000
