# MI VIS - Visual CMA

### Zadání 

Na datech s předchozího cvičení odskoušejte jednoducho verzi CMA. K výpočtu hodnot vzájemné informace můžete použít balíček Bioc2cor v R (dostupný v repozitáři CRAN).

Jako hodnocený úkol vytvořte program na bázi Pymolu, který pro libovolný protein v PDB, ke kterému je k dispozici vícečetné zarovnání (MSA) spočítá pozice se zajímavými hodnotami MI, tyto zobrazí a vytvoří pro ně tabulku s uvedením spočítané hodnoty a vzdálenosti příslušných aminokyselin v PDB modelu.

### Konfigurace

Konfigurace je obsažena v horní části programu

```python
FASTA_FILE    = "a.fasta"  # soubor s MSA sekvencemi nejvyrchnější musí být v PDB a komentář musí mít ve formátu >PDBID_CHAIN
TH_MI         = 0          # Práh pro MI která bude akceptována
TH_LEN        = 4          # Práh pro minimální délku amino kyselin sekvence
LEN_MAX       = 10         # Práh pro maximální délku amino kyselin sekvence
SHOW_ID       = 0          # id sekvence, která se zobrazí v pymol
SHOW_FIRST    = 10         # kolik sekvencí zobrazit v tabulce
MAX_SAME      = 0.95       # kolik % společné mohou mít jednotlivé sekvence <0, 1>
MI_F          = "mi"       # funkce pro spočítání vzájemné informace dostupné jsou mi a mip, mip je mi doplněna o normalizaci odstranujici evolucni informaci
```

### Spuštení

Program se spuští buď pomocí teminálu jako python `mi_vis.py` nebo přímo v pymol pomocí `spustit skript`. Pokud je spuštěn v pymol zobrazí i vybraný pár z konfigurace `SHOW_ID`, pokud je v terminálu pouze vypíše tabulku a uploží obrázky z pymol jednotlivých párů a skončí.

#### Závislosti

* Python
    * Pymol
    * numpy
    * matplotlib

* R
    * Bioc2cor

> `install.packages('Bios2cor', repos='http://cran.us.r-project.org')`

#### Popis tabulky

```txt
ID      START 1 START 2 SIZE    AVG MI  DISTANCE
0       8       13      4       0.0     6.82239294052124
1       6       17      5       0.0     16.43321990966797
2       10      17      5       0.0     11.177395820617676
3       4       18      3       0.0     19.320528030395508
4       11      18      5       0.0     11.021097183227539
5       5       19      4       0.0     19.609050750732422
6       12      19      5       0.0     11.004870414733887
7       11      20      5       0.0     13.71716594696045
8       10      21      5       0.0     16.2342472076416
9       12      21      5       0.0     13.56906509399414
```

tabulka se skládá ze 6 sloupců pricemz kazdy je odeleny tabulatorem. První sloupec obsahuje ID, které je vnitřním identifikátorem páru a odpovídá ID u obrázku. START1 je začátkem první párové sekvence a start2 je začátkem druhé párové sekvence v sekvenci protejnu. Size je délka podsekvencí a AVG MI je průměrná hodnota MI v této sekvenci. DISTANCE je vzálenost centeroidů těchto dvou sekvencí v protejnou z PDB.

### Způsob řešení

1. program vytvoří program v R, který použije knihovnu Bioc2cor a vypočítá korelační matici nad MSA sekvencemi a uloží ji do csv souboru
2. program v R je spuštěn pomocí os.system
3. Hlavní program přečte korelační matici z csv souboru
4. Hlavní program najde hodnoty v kovarienční matici které jsou meší než požadovaná hodnota z konfigurace
5. Program skouší následně tyto sekvence giagonálně rozšiřovat tak aby dodrřel podmínku na pořadované MI a pořadovanou Délku
6. Program vykreslí kovarienční matici s nalezenýmy sekvencemi
7. Pomocí vnitčních funkcí pymol obarví sekvence a změří mezi nimy vzdálenost pomocí mesure
8. Kařkou vyobrazenou dvojici uloží do PNG souboru
9. zobrazí vybranou dvojici a vypíše tabulku

