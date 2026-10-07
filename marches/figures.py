"""
Tracé des expériences avec matplotlib (facultatif).

    python -m marches.figures                       # liste les expériences
    python -m marches.figures grille_additive       # trace et affiche
    python -m marches.figures territoires_trois --sortie territoires.png

La simulation elle-même n'a besoin d'aucune bibliothèque ; seul ce module
demande matplotlib.
"""

import argparse

from .experiences import EXPERIENCES, executer, normaliser


COLONIES = ["Oranges", "Greens", "Purples"]


def tracer(sortie, axe=None):
    """Trace la sortie de `executer`, `simuler`, `simuler_1d` ou `balayer`."""
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    sortie = normaliser(sortie)

    if axe is None:
        _, axe = plt.subplots(figsize=(7, 7))

    if sortie["nature"] == "courbe":
        for serie in sortie["series"]:
            axe.errorbar(
                serie["valeurs"], serie["moyennes"],
                yerr=[e / serie["repetitions"] ** 0.5 for e in serie["ecarts_types"]],
                marker="o", capsize=2, label=serie["nom"],
            )

        axe.set_xlabel(sortie["variable"])
        axe.set_ylabel(sortie["libelle_mesure"])
        axe.legend()

    elif sortie["nature"] == "ligne":
        positions = sortie["positions"]
        points = list(zip(range(len(positions)), positions))

        _degrade(axe, LineCollection, points, "turbo_r")
        axe.set_xlabel("Temps (pas)")
        axe.set_ylabel("Position")
        axe.autoscale()

    else:
        taille = sortie["parametres"]["taille"]
        colonies = sortie["parametres"]["populations"] > 1

        for marche in sortie["marches"]:
            points = list(zip(marche["x"], marche["y"]))
            palette = COLONIES[marche["population"] % 3] if colonies else "turbo_r"

            _degrade(axe, LineCollection, points, palette)

        for x, y in sortie["departs"]:
            axe.plot(x, y, "o", color="black", markersize=6)

        axe.set_xlim(0, taille - 1)
        axe.set_ylim(0, taille - 1)
        axe.set_aspect("equal")

    axe.set_title(sortie.get("titre", ""))

    return axe


def _degrade(axe, LineCollection, points, palette) -> None:
    """Ligne brisée dont la couleur suit le temps."""
    if len(points) < 2:
        return

    segments = [points[i:i + 2] for i in range(len(points) - 1)]
    ligne = LineCollection(segments, cmap=palette, linewidths=0.8)

    # On évite les tons les plus clairs de la palette, peu lisibles.
    ligne.set_array([0.25 + 0.75 * i / len(segments) for i in range(len(segments))])
    ligne.set_clim(0, 1)
    axe.add_collection(ligne)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Trace une expérience du rapport."
    )
    parser.add_argument("experience", nargs="?")
    parser.add_argument("--graine", type=int)
    parser.add_argument("--sortie", help="fichier image ; sinon affichage à l'écran")
    args = parser.parse_args()

    if args.experience is None:
        for identifiant, experience in EXPERIENCES.items():
            print(f"{identifiant:<30}{experience['titre']}")
        return

    import matplotlib

    if args.sortie:
        matplotlib.use("Agg")

    import matplotlib.pyplot as plt

    tracer(executer(args.experience, graine=args.graine))

    if args.sortie:
        plt.savefig(args.sortie, dpi=150, bbox_inches="tight")
        print(f"Figure écrite dans {args.sortie}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
