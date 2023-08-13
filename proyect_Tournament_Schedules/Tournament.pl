% alan almora 189727
% Emiliano Mayen 192270

% -------------- Initialization populate -------------------
% Funtions:
% populate(i,i,i,o)
% setPopulation(i)
% games(i)
% sig(i,i,i,o,o)
% randomAssignment(i,o)
% approx(i,o)
% gameToAdd(i,i,o)
% available(i,i,o)
% count(i,i,o)


populate:-
    n(N),
    populate(0,N,[],Resp),
    setPopulation(Resp),
    !.
populate(N,N,Pop,Pop):-!.
populate(I,N,Temp,Pop):-
    I=<N,
    nTeams(Nteams),
    approx(Nteams,New),
    append([New],Temp,Pop2),
    I2 is I+1,
    populate(I2,N,Pop2,Pop),
    !.


setPopulation(Resp):-
    retractall(population),
    asserta(population(Resp):-!),
    !.

% Generate all the possibles combinations of matches
games(N):-
    games(1,1,N,Resp),
    write(Resp),
    !.
games(N,Res):-
    games(1,1,N,Res),
    !.
games(N,B,N,[[N,B]]):-
    B=:=(N-1),
    !.
games(A,A,N,Res):-
    sig(A,A,N,A2,B2),
    games(A2,B2,N,Res),
    !.
games(A,B,N,Res):-
    sig(A,B,N,A2,B2),
    games(A2,B2,N,Temp),
    append(Temp,[[A,B]],Res),
    !.
% Auxiliary function, Generate next values
sig(A,N,N,SigA,SigB):-
    SigA is (A+1),
    SigB is 1,
    !.
sig(A,B,_,A,SigB):-
    SigB is (B+1),
    !.

randomAssignment(N,Resp):-
    games(N,Temp),
    random_permutation(Temp,Resp),
    !.

% Generates a random games. Add games in a journee (if the teams are availables) 
approx(N,Resp):-
    games(N,Remaining),
    K is N/2,
    approx(Remaining,[],K,[],Resp),
    !.
approx(Remaining,_,_,Res,Out):-
    length(Remaining,Lr),
    Lr=:=0,
    Out=Res,
    !.
approx(Remaining,Week,K,Res,Out):-
    length(Week,L),
    L=:=K,
    approx(Remaining,[],K,Res,Out),
    !.
approx(Remaining,Week,K,Res,Out):-
    length(Remaining,Lr),
    Lr>0,
    gameToAdd(Remaining,Week,Game),
    % Add to Res
    append(Res,[Game],Res2),
    % Add to Week
    append([Game],Week,Week2),
    % Remove from remaining
    subtract(Remaining,[Game],Remaining2),
    % Recursive call
    approx(Remaining2,Week2,K,Res2,Out),
    !.

% game to add in a journee  
gameToAdd(Remaining,Week,Game):-
    available(Remaining,Week,Available),
    length(Available,La),
    La>0,
    random_member(Game,Available),
    !.
gameToAdd(Remaining,_,Game):-
    random_member(Game,Remaining),
    !.
    
% count games by journee
gameTeamsInWeek(TeamsInWeek,[A|[B]]):-
    count(TeamsInWeek,A,Ac),
    Ac=:=0,
    count(TeamsInWeek,B,Bc),
    Bc=:=0,
    !.
available(Remaining,Week,Res):-
    flatten(Week,TeamsInWeek),
    include(gameTeamsInWeek(TeamsInWeek),Remaining,Res),
    !.

% Auxiliary function, count the number of occurances in a list.
count([],_,0).
count([H|T],H,NewCount):-
    count(T,H,OldCount),
    NewCount is OldCount +1,
    !.
count([H|T] , H2,Count):-
 dif(H,H2),
 count(T,H2,Count),
 !.

% -------------------- Fitnsess --------------------
% Funcitons:
% fitness(i,o).



% Determines the fitness of an individual.
fitness(Individual,F):-
    nTeams(N),
    conflictos(Individual,N,Conf),
    length(Individual,MaxConf),
    onlyValidScore(OnlyValid),
    (
        OnlyValid-> Alt is 1; alterningScore(Individual,N,Alt)
    ),
    Temp is (0.95*(1-(Conf/MaxConf)))+(Alt*0.05),
    (
        Conf=:=0 -> Bonus is 1.0;
        Bonus is 0.0
    ),
    F is Temp+Bonus,
    !.

conflictos(Assignment,N,Out):-
    conflictos(Assignment,0,N,0,WeekConf),
    list_to_set(Assignment,Set),
    length(Set,S),
    A is (N*(N-1)),
    Out is WeekConf+((A-S)/A),
    !.
conflictos(_,I,N,Count,Out):-
    I>=(N*(N-1)),
    Out=Count,
    !.
conflictos(Assignment,I,N,Count,Out):-
    End is I+(N/2)+1,
    slice(Assignment,I,End,Slice),
    flatten(Slice,F),
    list_to_set(F,Set),
    length(Set,S2),
    Count2 is Count+(N-S2)/(N/2),
    I2 is (End-1),
    conflictos(Assignment,I2,N,Count2,Out),
    !.

% Rewards an assignment for having teams alternate between
% home and visiting games.
% alterningScore(i,i,o).
alterningScore(Assignment,N,Out):-
    alterningScore(Assignment,1,N,Sum),
    Out is Sum/N,
    !.

alterningScore(_,Team,N,Sum):-
    Team>N,
    Sum is 0,
    !.
alterningScore(Assignment,Team,N,Sum):-
    Team2 is (Team+1),
    alterningScore(Assignment,Team2,N,Temp),
    teamGames(Assignment,Team,[[A|[_]]|Rest]),
    (
        A=:=Team -> Home is 1; Home is 0
    ),
    length(Rest,L),
    Ngames is (L-1),
    alternar(Rest,Team,Home,Alt),
    (
        Ngames=:=0 -> Sum is Temp;
        Sum is Temp+((Alt-1)/Ngames)
    ),
    !.

teamInGame(A,[A|[_]]):-!.
teamInGame(B,[_|[B]]):-!.
teamGames(Assignment,Team,Res):-
    include(teamInGame(Team),Assignment,Res),
    !.

alternar([],_,_,0):-!.
alternar([[_|[B]]|Rest],B,1,Alt):-
    alternar(Rest,B,0,Temp),
    Alt is Temp+1,
    !.
alternar([[A|[_]]|Rest],A,0,Alt):-
    alternar(Rest,A,1,Temp),
    Alt is Temp+1,
    !.
alternar([[A|[_]]|Rest],A,1,Alt):-
    alternar(Rest,A,1,Alt),
    !.
alternar([[_|[B]]|Rest],B,0,Alt):-
    alternar(Rest,B,0,Alt),
    !.



% ----------------------- main -----------------------

:- dynamic population /1. 
:- dynamic listOfGames /1.
:- dynamic coldness/1.

library(lists).
library(apply).
library(random).
library(pairs).


% Parameters 
n(30). % Size of population
nTeams(16). % Number of teams in tournament
propElite(0.05). % Best % of population considered elite
mutationRate(0.2). % Rate at which to randomly apply mutation
coldness(0.5). % Selection pressure parameter used in tournament selection
coldnessInc(0.005). % Amount with which to increment coldness each generation.
maxColdness(0.85). % The max value coldness can take.
% Used in stopping criteria. Max number of iterations allowable without improvement.
maxWithoutImprovement(10). 
onlyValidScore(true).


% Derived parameters
% Number of individuals of elite as a function of n and propElite.
% elites(o).
elites(Res):-
    n(N),
    propElite(Prop),
    P is N*Prop,
    Res is floor(P),
    !.
% Number of non-elite individuals.
% noElites(o).
noElites(Res):-
    elites(E),
    n(N),
    Res is (N-E),
    !.

% Main predicate. Initializes population, evolves them and outputs solutions.
evoluciona:-
    % Write parameters
    parameters,

    % Set list of games to mutate faster
    nTeams(Nteams),
    games(Nteams,G),
    setGames(G),

    % Initialize population
    populate,
    write("Initialized population with approx solutions"),nl,nl,

    % evoluciona
    maxWithoutImprovement(M),
    evoluciona(1,0,-100,M),

    % Output solutions
    setOfValid(Out),
    length(Out,L),
    L>0,
    nl,write("Found "),write(L),write(" unique valid solution(s)."),nl,
    write("Writing them to output.txt"),nl,
    outputToFile("output.txt"),
    !.

% If failed to find solution, try again.
evoluciona:-
    write("Trying again..."),nl,nl,
    evoluciona,

    !.

% Main evoluciona loop. Loops while TimeSince<M (maxWithoutImprovement).
% evoluciona(i,i,i,i).
evoluciona(_,M,_,M):-!.
evoluciona(Gen,TimeSince,LastMax,M):-
    % Set aside elite
    elite(Elite),

    % Apply selection, crossover and mutation
    tournamentSelection,
    crossoverPopulation,
    mutarPoblacion,

    % Append elite to population
    population(Pop),
    append(Elite,Pop,New),
    setPopulation(New),

    % Calculate fitness of individuals
    fitnesses(New,Fs),

    % Determine Time SInce Improvement
    max_list(Fs,CurrentMax),
    newTimeSinceImprovement(TimeSince,LastMax,CurrentMax,TimeSince2,Max2),

    % Write relevant generation statistics
    writeStatistics(Gen,Fs,TimeSince2),

    % Update selection pressure
    updateColdness,

    % Repeat
    Gen2 is (Gen+1),
    evoluciona(Gen2,TimeSince2,Max2,M),
    !.

updateColdness:-
    coldness(Old),
    coldnessInc(Inc),
    maxColdness(M),
    New is (Old+Inc),
    Actual is min(New,M),
    setColdnesss(Actual),
    !.


% newTimeSinceImprovement(i,i,i,o,o).
newTimeSinceImprovement(_,LastMax,Max,NewTime,NewMax):-
    Max>LastMax,
    NewTime is 0,
    NewMax is Max,
    !.
newTimeSinceImprovement(LastTime,LastMax,_,NewTime,NewMax):-
    NewTime is (LastTime+1),
    NewMax is LastMax,
    !.


% elite(o).
elite(Res):-
    elites(N),
    population(Pop),
    sortByFitness(Pop,Sorted),
    firstN(Sorted,N,Res),
    !.

% Tournament Selection.
% Randomly draws 2 individuals and a random number.
% If random number < coldness (selection pressure),
% keep the best individual (according to fitness),
% else keep the worst.
% Repeat NoElites times.
% tournamentSelection.
tournamentSelection:-
    noElites(N),
    tournamentSelection(0,N,[],New),
    setPopulation(New),
    !.
% tournamentSelection(i,i,i,o).
tournamentSelection(N,N,Temp,Temp):-!.
tournamentSelection(I,N,Temp,New):-
    population(Pop),
    random_member(A,Pop),
    random_member(B,Pop),
    mejoryPeor(A,B,Mejor,Peor),
    random(0.0,1.0,R),
    toAdd(Mejor,Peor,R,ToAdd),
    append([ToAdd],Temp,Temp2),
    I2 is (I+1),
    tournamentSelection(I2,N,Temp2,New),
    !.

% Used in tournament selection. Determines if the worse or
% better individual is added as a function of temperature (selection pressure),
% and the drawn random number.
% toAdd(i,i,i,o).
toAdd(Mejor,_,R,ToAdd):-
    coldness(T),
    R<T,
    ToAdd=Mejor,
    !.
toAdd(_,Peor,_,Peor):-!.

% Mutate the population.

mutarPoblacion:-
    population(Pop),
    mutarPoblacion(Pop,[],New),
    setPopulation(New),
    !.
% mutarPoblacion(i,i,o).
mutarPoblacion([],Temp,Temp):-!.
mutarPoblacion([I|Rest],Temp,New):-
    random(0.0,1.0,R),
    mutationRate(Mr),
    R<Mr,
    mutar(I,Mutated),
    append(Temp,[Mutated],Temp2),
    mutarPoblacion(Rest,Temp2,New),
    !.
mutarPoblacion([I|Rest],Temp,New):-
    append(Temp,[I],Temp2),
    mutarPoblacion(Rest,Temp2,New),
    !.

% Mutation operator.
% Replace an individual's  game at a random index with a random game
% taken from games.
% mutar(i,o).
mutar(Individual,Mutated):-
    listOfGames(Gs),
    random_member(Game,Gs),
    length(Individual,U),
    random(0,U,Index),
    remplaza(Individual,Index,Game,Mutated),
    !.


% Crossover the population.
% Crossover two random memebers of the population and add their offspring
% to the population. Repeat NoElites/2 times.
crossoverPopulation:-
    noElites(N),
    Half is N/2,
    crossoverPopulation(0,Half,[],New),
    setPopulation(New),
    !.
% crossoverPopulation(i,i,i,o).
crossoverPopulation(I,N,Temp,Temp):-
    I>=N,
    !.
crossoverPopulation(I,N,Temp,New):-
    population(Pop),
    random_member(A,Pop),
    random_member(B,Pop),
    crossover(A,B,Off1,Off2),
    append([Off1,Off2],Temp,Temp2),
    I2 is (I+1),
    crossoverPopulation(I2,N,Temp2,New),
    !.
    
% Linear 1 point crossover
% crossover(i,i,i,o,o):
crossover(A,B,Off1,Off2):-
    length(A,U),
    Upper is (U+1),
    random(1,Upper,Index),
    crossover(A,B,Index,Off1,Off2),
    !.
crossover(L1, L2, Index, H1, H2):-
    IM1 is Index-1,
    length(L1, Length1),
    slice(L1, 0, Index, L1P1),
    slice(L2, 0, Index, L2P1),
    Length is Length1+1,
    slice(L1, IM1, Length, L1P2),
    slice(L2, IM1, Length, L2P2),
    append(L1P1, L2P2, H1),
    append(L2P1, L1P2, H2).


% ------------ Delate -----------
% Negative of finess. Used for sorting.
% negFiness(i,o).
negFitness(Individual,F):- %For sorting
    fitness(Individual,Temp),
    F is Temp*(-1),
    !.

% Determines the better and worse from individuals A and B.
% mejoryPeor(i,i,o,o).
mejoryPeor(A,B,Mejor,Peor):-
    fitness(A,Fa),
    fitness(A,Fb),
    Fa>Fb,
    Mejor is A,
    Peor is B,
    !.
mejoryPeor(A,B,B,A):-!.

% Sort by fitness in descending order.
% sortByFitness(i,o).
sortByFitness(L,Sorted):-
    population(L),
    map_list_to_pairs(negFitness,L,Pairs),
    keysort(Pairs,Temp),
    pairs_values(Temp, Sorted),
    !.

% Detrmines de fitness of a list of individuals.
% fitnesses(i,o).
fitnesses(Pop,Out):-
    fitnesses(Pop,[],Out),
    !.
fitnesses([],Temp,Temp):-!.
fitnesses([I|Rest],Temp,Out):-
    fitness(I,F),
    append([F],Temp,Temp2),
    fitnesses(Rest,Temp2,Out),
    !.

% Determines if an assignment is valid (conflictos is 0).
% valid(i).
valid(Assignment):-
    nTeams(N),
    conflictos(Assignment,N,Out),
    Out=:=0,
    !.

% Determines the set of individuals from population that are valid.
% setOfValid(o).
setOfValid(Res):-
    population(Pop),
    include(valid,Pop,List),
    list_to_set(List,Res),
    !.


% Slice a list [From,To).
% slice(i,i,i,o).
slice(L, From, To, R):-
  length(LFrom, From),
  length([_|LTo], To),
  append(LTo, _, L),
  append(LFrom, R, LTo).

% Determines first N elements of a list.
% fitstN(i,i,o).
firstN(L,N,Res):-
    End is (N+1),
    slice(L,0,End,Res),
    !.

% Determines last N elements of a list.
% last(i,i,o).
lastN(L,N,Res):-
    length(L,End),
    ActualEnd is (End+1),
    S is (End-N),
    ActualStart is max(S,0),
    slice(L,ActualStart,ActualEnd,Res),
    !.

% Replaces the element at an index with a new element.
% remplaza(i,i,i,o).
remplaza([_|T], 0, X, [X|T]):-!.
remplaza([H|T], I, X, [H|R]):- I > -1, NI is I-1, remplaza(T, NI, X, R), !.
remplaza(L, _, _, L).



% Sets games as L in dynamic database.
% setGames(i).
setGames(L):-
    retractall(listOfGames),
    asserta(listOfGames(L):-!),
    !.

% Sets coldness as C in dynamic database.
% setColdnesss(i).
setColdnesss(C):-
    retractall(coldness),
    asserta(coldness(C):-!),
    !.

% Writes Parameters
parameters:-
    write("\nShowing parameters"),nl,nl,
    n(N),
    write("Population size: "),write(N),nl,
    nTeams(Nteams),
    write("Number of teams: "),write(Nteams),nl,
    mutationRate(Mr),
    write("Mutation rate: "),write(Mr),nl,
    propElite(P),
    elites(E),
    write("Elite proportion: "),write(P),write(" ("),write(E),write(")"),nl,
    !.

% Writes relevant statistics.
% writeStatistics(i,i,i).
writeStatistics(Gen,Fs,TimeSince):-
    write("Gen: "),write(Gen),
    min_list(Fs,Min),
    max_list(Fs,Max),
    average(Fs,Mean),
    setOfValid(Valid),
    coldness(C),
    length(Valid,V),
    write("   Min: "),format("~4f", [Min]),
    write("   Mean: "),format("~4f", [Mean]),
    write("   Max: "),format("~4f", [Max]),
    write("   # valid: "),write(V),
    write("   Time Since Improvement: "),write(TimeSince),
    nl,  
    !.


% Determines the average of a list.
% average(i,o).
average(List,Average):- 
    sumlist(List,Sum),
    length(List,Length),
    Length>0, 
    Average is Sum/Length.



outputToFile(Filename):-
    setOfValid(List),
    open(Filename, write, File),
    writeList(File, List),
    close(File).

writeList(_File, []) :- !.
writeList(File, [Head|Tail]) :-
    write(File, Head),
    fitness(Head,F),
    write(File,","),write(File,F),
    write(File, '\n'),
    writeList(File, Tail).