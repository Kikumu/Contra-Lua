# Contra-Lua.
--Genetic algorithm implemented using lua and bizhawk to play classic contra retro game--

The goal of this project is to test my understanding of the NEAT algorithm. I got very heavy inspiration from sethblings MAR/IO video and code on NEAT algorithms alongside a
number of research papers which I am currently looking into and will list/credit in this read me.

**Please note: This project is still undergoing and I intend to complete it as soon as i get some free time!**

# Roadmap part 1.

1) Get to grips with Lua programming language **DONE**
2) Understand how to use TAS and Bizhawk to extract game information **DONE**
3) Extract game data including player health, locations, enemy status **DONE**
4) In conjunction with looking at research papers and seth blings code, rework and understand the code to fit the solution of my current problem **DONE**


# Roadmap part 2(Mapping Damage states/ Bug squashing in progress).

1) Create connection gene structures **DONE**
2) Create node gene structure **DONE**
3) Create sample network **DONE**
4) Create mutation functions **DONE** 
5) Test mutation functions **DONE**
6) Test network with mutation functions **DONE**
7) Create crossover functions **DONE**
8) Test crossover functions **DONE**
9) Test crossover, mutation and networ functions **DONE**
10) Create species function **DONE**
11) Test species function**DONE**
12) Create fitness calculation functions**DONE**
13) Test fitness calculation functions**DONE**
14) Run and evaluate sample generation**DONE**
15) Create stale species function(remove species which arent improving) **DONE**
16) Create function which allows only species which are doing well to breed more and not doing well to breed less.(what i have atm is that, all species are allowed to bread equivalently no matter the fitness this not very nice)**DONE**
15) Official test run in game **NOT DONE**


# Resources.
1) http://nn.cs.utexas.edu/downloads/papers/stanley.ec02.pdf
2) https://www.researchgate.net/publication/223999345_Tree-Based_Indirect_Encodings_for_Evolutionary_Development_of_Neural_Networks
3) https://www.youtube.com/watch?v=qv6UVOQ0F44 (Seth's marl/io)
4) http://www.cs.ucf.edu/~kstanley/neat.html

