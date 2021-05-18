## Istanza d'esempio

| t | 1  | 2 | 3  | 4  | 5  |
|---|----|---|----|----|----|
| s | 50 | 0 | 30  | 50 | 0  |
| B | 50 | 0 | 30  | 20 | 0  |
| x |    |   | 50 |    |    |
| y |    |   |    |    | 50 |

$\Delta = 2$

$t = T - \Delta = 5 - 2 = 3$

$c = min\{s(t+\Delta-1); x(t); B(t-1) \} = min\{s(4);x(3); \textbf{B(2)}\} = min\{50;50;\textbf{0} \}$ = 0

...

$x(t-1) = x(2) + c = 0$

...

$t = t-1 = 2$, ma ho che $x(2) = 0$, dunque quando vado a calcolare il valore di c, avrò sempre 0 come valore, il che si ripercuote sul valore degli altri valori di $x$ che rimane a $0$.



