# Validation camerounaise de la méthode de Vannes

Ce dépôt conserve la transcription de la méthode point par point de Nathalie
Vannes (2011) et développe une méthode dérivée physiquement cohérente pour la
correction classique de terrain : la **condensation résiduelle centrée sur la
station** (CRCS).

Sur les 1 991 stations, dans des conditions strictement identiques, la CRCS
présente par rapport aux prismes exacts de Nagy un biais de -0,002618 mGal, une
MAE de 0,002619 mGal, une RMSE de 0,006768 mGal et une corrélation de 0,999999.
La formulation originale, utilisée directement comme correction de terrain,
présente une RMSE de 5,098439 mGal et 604 valeurs négatives.

Ce résultat rend la méthode corrigée **crédible comme équivalent numérique de
Nagy dans ce protocole contrôlé**, mais pas encore comme validation indépendante
définitive : les cinq premiers kilomètres sont régularisés par les prismes de
Nagy et la grille fournie est interpolée, finie et sans paramètres de krigeage
reproductibles.

Les équations, tests, limites et conditions nécessaires à une publication sont
documentés dans le dossier `docs/`.

## Données, droits et diffusion

Le logiciel est placé sous licence MIT. Cette licence ne couvre pas les données.
Cameroon Mining Group S.A. a autorisé la publication et l'archivage scientifique
permanent des coordonnées pseudonymisées, altitudes de station, grilles de
terrain, corrections calculées, figures et statistiques utilisées dans cette
étude (autorisation 00352/10/09/26/CMG). La propriété des données sources reste
à Cameroon Mining Group S.A. Les mesures gravimétriques brutes, noms de clients
ou de prospects, interprétations commerciales et métadonnées opérationnelles
sont exclus. Voir `DATA_USE_NOTICE.md`.

Le classeur public des stations ne contient donc que la longitude, la latitude
et l'altitude. La colonne d'anomalie d'air libre de la version de travail interne
est volontairement exclue de la version `v1.0.0`.

## Auteur unique

Ariel Bounou Guimking — Faculty of Physics and Technology, Novosibirsk State
Technical University; ORCID : [0009-0003-1847-2436](https://orcid.org/0009-0003-1847-2436).

