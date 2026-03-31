# Exemple Flowshop (atelier cyclique)

Periodic ﬂowshop.
- Parts are carried on pallets. When the tasks on a part are ﬁnished the pallet start another cycle with another part of the same class. Each part visit machines in sequence never coming back to the same machine.
- To make a task we can have many machines (called a class). The machines visit the parts in sequence never coming back to the same part.
- We deﬁne the ﬂowshop by a matrix describing the resources used and the processing times.
- Each line of the matrix is associated to a machine class.
- Each columns to a part class.
- The entries of the matrix are the processing times.
- If a part class does not need a machine class the corresponding entry is −∞

## Données : temps de traitement `PT`

Machines × pièces en `(max,+)`. Une case vide correspond à `mp0` (`ε`).

```julia
using MaxPlus

const PT = MP.([
    2 3.9 0.95 1.1 0.7 1.4
    mp0 mp0 2 1.2 mp0 1.7
    3.7 mp0 2.2 mp0 6.4 mp0
    mp0 mp0 2 mp0 1 1
    1.7 3.1 3 mp0 1.3 mp0
    0.5 3.2 4.3 1.9 1.6 0.4
    1 1 1 1 1 1
    1.5 1.5 1.5 1.2 1.2 1.2
])

nmach, npiece = size(PT)  # (8, 6)

# Nombre de machines / palettes par « classe » (ici 1 partout — cas périodique 1)
nm = ones(Int, nmach)
np = ones(Int, npiece)
```

```julia
T, N = flowshop_graph(PT, Float64.(nm), Float64.(np))
r = mpeigen(T, N) # Semi Howard
λ = plustimes(r.eigenvalues[1])  # cadence-type (semi-Howard)
```

## Système linéaire et simulation en boucle

```julia
s = flowshop(PT)
nt = 100
u = ones(Float64, nmach + npiece, nt)   # colonnes = temps, lignes = entrées

chi, y = flowshop_simu(s, nm, np, u)
# y : Float64, taille (nt, nmach + npiece), dérive cyclique déjà soustraite
```

**ScicosLab (équivalences)** : `fbm`, `fbp` / `shift`, composition `|`, feedback `/`, puis `explicit` et `simul` sont encapsulés dans la fonction Julia `flowshop_simu`.

## Tracés (Plots.jl)

Les scripts ci-dessous reproduisent les courbes affichées plus bas.

```julia
using Plots
Plots.default(size = (760, 400), legend = :outerright)

p = plot(title = "Transitoire, m = p = 1", xlabel = "k", ylabel = "y")
for j in 1:min(14, size(y, 2))
    plot!(p, 1:size(y, 1), y[:, j]; label = "sortie $j")
end
savefig(p, "flowshop_transient_m1.svg")

# Cas « buffers » 3 (refaire la simulation avec d’autres nm, np)
nm3 = fill(3, nmach)
np3 = fill(3, npiece)
_, y3 = flowshop_simu(s, nm3, np3, u)
p2 = plot(title = "Transitoire, m = p = 3", xlabel = "k", ylabel = "y")
for j in 1:min(14, size(y3, 2))
    plot!(p2, 1:size(y3, 1), y3[:, j]; label = "sortie $j")
end
savefig(p2, "flowshop_transient_m3.svg")
```

### Résultat — périodicité 1 (`m = p = 1`)

![Transitoire flowshop, m=p=1](assets/flowshop_transient_m1.svg)

### Résultat — buffers 3 (`m = p = 3`)

Dans le script Scilab d’origine, le second tracé réutilisait parfois `y` sans relancer la simulation ; ici on **recalcule** `flowshop_simu` avec `nm = 3`, `np = 3`.

![Transitoire flowshop, m=p=3](assets/flowshop_transient_m3.svg)
