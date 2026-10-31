# Validation de l'adaptation de Vannes centrée sur la station

Ce calcul contrôlé utilise les 1,991 stations et la totalité de la
grille historique informationnellement équivalente de 200 x 200 nœuds. Toutes
les méthodes reçoivent les mêmes cellules, la même densité, la même altitude de
station et le même traitement des altitudes négatives.

| Variante comparée à Nagy | Biais (mGal) | MAE (mGal) | RMSE (mGal) | P95 absolu (mGal) | Maximum absolu (mGal) | Pearson r | Valeurs négatives |
|---|---:|---:|---:|---:|---:|---:|---:|
| Vannes résiduelle originale | -1.866119 | 2.323545 | 5.098439 | 9.448417 | 73.950194 | -0.274040 | 604 |
| Adaptation centrée sur la station | -0.002618 | 0.002619 | 0.006768 | 0.011113 | 0.127825 | 0.999999 | 0 |

Pour la méthode adaptée, 93.97 % des stations diffèrent de Nagy
de moins de 0,01 mGal, 99.75 % de moins de 0,05 mGal et
99.95 % de moins de 0,10 mGal.

L'adaptation n'est ni une régression ni un réétalonnage empirique. Elle déplace
la surface de condensation de référence du géoïde vers la sphère passant par la
station, afin que la masse résiduelle et la masse de référence représentent la
même colonne. Les cellules proches sont calculées par les prismes exacts de Nagy
afin de régulariser la cellule contenant la station.

Ce résultat valide l'orientation mathématique de la correction sur la grille
historique fournie. Les valeurs destinées à une publication nécessitent encore
un nouveau calcul de la chaîne MNT hybride/FABDEM et de sa convergence radiale,
avec la colonne `vannes_station_centred_mgal` activée.
