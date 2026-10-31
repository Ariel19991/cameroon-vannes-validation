# Analyse de sensibilité

L’analyse utilise 80 stations réparties régulièrement dans l’ordre du fichier.
À pas de grille constant (40 000 cellules), les ordres de quadrature 4, 6 et 8
donnent la même RMSE à mieux que (10^{-9}) mGal : l’intégration radiale est
numériquement convergée dès l’ordre 4 pour ce jeu de données.

| Pas/stride | Champ proche | Ordre | RMSE à Nagy (mGal) | P95 absolu (mGal) | Maximum (mGal) |
|---:|---:|---:|---:|---:|---:|
| 3 | 0 km | 4 | 0,012425 | 0,013783 | 0,074269 |
| 3 | 2 km | 4 | 0,012425 | 0,013783 | 0,074269 |
| 3 | 5 km | 4 | 0,006764 | 0,006587 | 0,039001 |
| 3 | 10 km | 4 | 0,001050 | 0,002264 | 0,004520 |
| 6 | 5 km | 4 | 0,014365 | 0,018995 | 0,084792 |
| 12 | 5 km | 4 | 0,007380 | 0,015954 | 0,021347 |

La baisse de l’écart lorsque le rayon Nagy augmente est attendue et ne constitue
pas une preuve indépendante : une fraction croissante du résultat devient
identique à la référence par construction. En revanche, la stabilité aux ordres
4–8 montre que l’erreur restante ne provient pas de la quadrature radiale. La
non-monotonie avec le stride confirme que la discrétisation et le MNT dominent
l’incertitude résiduelle.

