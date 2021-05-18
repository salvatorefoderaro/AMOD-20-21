##  Istanza d'esempio
  
  
| t | 1  | 2 | 3  | 4  | 5  |
|---|----|---|----|----|----|
| s | 50 | 0 | 30  | 50 | 0  |
| B | 50 | 0 | 30  | 20 | 0  |
| x |    |   | 50 |    |    |
| y |    |   |    |    | 50 |
  
<img src="https://latex.codecogs.com/gif.latex?&#x5C;Delta%20=%202"/>
  
<img src="https://latex.codecogs.com/gif.latex?t%20=%20T%20-%20&#x5C;Delta%20=%205%20-%202%20=%203"/>
  
<img src="https://latex.codecogs.com/gif.latex?c%20=%20min&#x5C;{s(t+&#x5C;Delta-1);%20x(t);%20B(t-1)%20&#x5C;}%20=%20min&#x5C;{s(4);x(3);%20&#x5C;textbf{B(2)}&#x5C;}%20=%20min&#x5C;{50;50;&#x5C;textbf{0}%20&#x5C;}"/> = 0
  
...
  
<img src="https://latex.codecogs.com/gif.latex?x(t-1)%20=%20x(2)%20+%20c%20=%200"/>
  
...
  
<img src="https://latex.codecogs.com/gif.latex?t%20=%20t-1%20=%202"/>, ma ho che <img src="https://latex.codecogs.com/gif.latex?x(2)%20=%200"/>, dunque quando vado a calcolare il valore di c, avrò sempre 0 come valore, il che si ripercuote sul valore degli altri valori di <img src="https://latex.codecogs.com/gif.latex?x"/> che rimane a <img src="https://latex.codecogs.com/gif.latex?0"/>.
  
  
  
  