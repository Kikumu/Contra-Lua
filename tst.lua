--[[
TO DO:
1) Take it out for a test spin to test out functions(mutations,crossovers,propagation,disjoint genes, matching genes, excess,outputs, inputs,network feedforward monitor) created --DONE
2) Test propagation --Done
3) Neurons holder for each network created --DONE
4) Find a way to keep track of the fitness of each genome --done
5) Create function to organise a group of genes into different speciesPool --in progress(speciation)
6) Fitness function of genomes
7) Rework how the network breeds --done
8) Fix how genes mutate --done
9) Integrate to game
10) crossover hard test......done
]]
local InitialPopulation1 = {}
local TestInputs = {}
TestInputs[1] = 1
TestInputs[2] = 3
TestInputs[3] = 4
TestInputs[4] = 6
TestInputs[5] = 7

local TestOutputs = {}
TestOutputs[1] = 0.8
TestOutputs[2] = 0.8


local TestOutputs1 = {}
TestOutputs1[1] = 10
TestOutputs1[2] = 10


--stores all species
speciesStore = {}
--print("Species first :"..#speciesStore)
--used to create species
--GEN --SPECIES--GENOME
--geneActivationChance = 0.3
genomeFlipChance = 0.45
genomeDecrementChance = 0.45
genomeMutationChance = 0.52
genomeActivationChance = 0.45
genomeLinkMutationChance = 0.35
genomeNodeMutationChance = 0.41
genomeStepValue = 0.55

--genomePointMutateChance = 0.4
selectionChance = 0.3 --in use
initialPopulationSize = 150 --in use
crossOverChance = 0.75
NumberOfGenerations = 100
--species classification variables
disjointmentConstant = 0.25 --c1
excessGenesConstant = 0.3 --c2
weightImportanceConstant = 0.25 --c3
speciesDistance = 0.70 --in use
MaxNodes = 100
MaxLinks = 30
--E is number of disjointed connection genes
--D is number of excess connection genes
--W is weight value
--N number of connection genes


local InitialPopulation = {}


function softMax(Outputs)
  --val = val/sum of exp(val) entire val output net including val itself
  end

--activation function
function sigmoid(x)
    return 1.0 / (1.0 + math.exp(-x))
end

function connectionGene()
local gene = {}
gene.input = 0
gene.out = 0
gene.weight = 0
gene.status = true
gene.innovation = 0 --ancestry monitor
return gene
end

--tested(needs reviewing) --TESTED AND REWORKED
function mutateConnectionGene(genome)
v1 = math.random(1,#genome.network)
v2 = math.random(1,#genome.network)
node1 = genome.network[v1] --1
node2 = genome.network[v2] --1
if genome.network[v1].inStatus == 0 and genome.network[v2].inStatus == 0 then
  print("neurons are both input neurons")
  return 0
end
if genome.network[v1].inStatus == 1 and genome.network[v2].inStatus == 1 then
  print("neurons are both output neurons")
  return 0
end
for key in genome.network[v1].geneInnovationIndex do
  if genome.network[v2].geneInnovationIndex[key]~=nil then
    print("Found matching gene")
    return  0
  end
end
if genome.network[v1].inStatus ~= genome.network[v2].inStatus then
  connectionGeneMutate1 = connectionGene()
  connectionGeneMutate1.innovation = genome.maxInnovation + 1
  genome.maxInnovation = connectionGeneMutate1.innovation
  v1 = false
  v2 = false
  if genome.network[v1].inStatus == 0 or genome.network[v2].inStatus == 2 then
    connectionGeneMutate1.in = v1
    v1 = true
  elseif genome.network[v2].inStatus == 0 or genome.network[v1].inStatus == 2 then
    v2 = true
    connectionGeneMutate1.in = v2
  end
  
  if (genome.network[v1].inStatus == 1 or genome.network[v2].inStatus == 2) and v1 == false then
    connectionGeneMutate1.out = v1
  elseif (genome.network[v2].inStatus == 1 or genome.network[v1].inStatus == 2) and v2 == false then
    connectionGeneMutate1.out = v2
  end
  
  genome.genes[connectionGeneMutate1.innovation] = connectionGeneMutate1
  genome.maxInnovation = connectionGeneMutate1.innovation
  end
end

--FULLY REWORKED
function mutateNodeGene(genome)
tstval = #genome.genes
maxInnovation = genome.maxInnovation
maxNeuron = genome.maxNeuron
maxNeuron = maxNeuron + 1
Neuron = newNeuron()
--print("Mutation max prev innovation value: "..maxInnovation)
genepos = math.random(1,#genome.genes)
genome.genes[genepos].status = false

connectionGeneMutate1 = connectionGene()
connectionGeneMutate1.input = genome.genes[genepos].input --attatch neuron
connectionGeneMutate1.innovation = maxInnovation + 1
connectionGeneMutate1.weight = 1
connectionGeneMutate1.out = maxNeuron


connectionGeneMutate2 = connectionGene()
connectionGeneMutate2.input = maxNeuron --attatch neuron
connectionGeneMutate2.innovation = connectionGeneMutate1.innovation + 1
connectionGeneMutate2.weight = genome.genes[genepos].weight
connectionGeneMutate2.out = genome.genes[genepos].out


Neuron.geneInnovationIndex[connectionGeneMutate1.innovation] = connectionGeneMutate1.innovation
Neuron.geneInnovationIndex[connectionGeneMutate2.innovation] = connectionGeneMutate2.innovation
genome.network[genome.genes[genepos].input].geneInnovationIndex[connectionGeneMutate1.innovation] = connectionGeneMutate1.innovation
genome.network[genome.genes[genepos].out].geneInnovationIndex[connectionGeneMutate2.innovation] = connectionGeneMutate2.innovation


genome.genes[connectionGeneMutate1.innovation] = connectionGeneMutate1
genome.genes[connectionGeneMutate2.innovation] = connectionGeneMutate2
genome.network[maxNeuron] = Neuron


genome.maxInnovation = connectionGeneMutate2.innovation
genome.maxNeuron = maxNeuron
end


function pointMutateGenome(genome)
  for i = 1, #genome.genes do
    if genome.Perturbance > math.random() then
      genome.genes[i].weight = genome.genes[i].weight + math.random() * genome.step*2.0 - genome.step
    else
      genome.genes[i].weight = math.random() * 4.0 - 2.0
    end
end
end

function matchingGenes(genome1,genome2)
  matches = 0
  sum_ = 0
  for innovation in genome1.genes do
    if genome2.genes[innovation]~=nil then
      if selectionChance > math.random() then
        matches = matches + 1
        sum_ = sum_ + math.abs(genome1.genes[i].weight - genome2.genes[i].weight)
      else
        matches = matches + 1
        sum_ = sum_ + math.abs(genome1.genes[i].weight - genome2.genes[i].weight)
        end
      end
    end
  average = sum_/matches
  return matches,average
end



function disjointGenes(genome1,genome2)
max = genome1.maxInnovation
max2 = genome2.maxInnovation
maximum = math.max(max,max2)
disjoints = 0

for innovation in genome1.genes do
  if genome2.genes[innovation]==nil and innovation <= genome1.maxInnovation then
    disjoints = disjoints + 1
  end
end
return disjoints
end

function excessGenes(genome1,genome2)
max = genome1.maxInnovation
excess = 0
  for innovation in genome2.genes do
    if innovation > max then
      excess = excess + 1
      end
  end
return excess
end


function crossover(genome1,genome2) --genome in the sense that you are passing to this function a set of genomes


if genome.mutationChance > math.random() then
--print("Point mutate attempt")
pointMutateGenome(genome)
--print("pintMutate")
end
if genome.linkMutationChance > math.random() then
--print("ConnectMutate Attempt")
mutateConnectionGene(genome)
end
if genome.nodeMutationChance > math.random() then
--print("NodeMutate Attempt")
mutateNodeGene(genome)
end

--print("Crossover genes number after mutations: "..#genome.genes)
return genome
end

function createNewGenome()
  genome = {} --holds the gene information network (general info)
  genome.genes = {} --weight info (connection genes info)
  genome.fitness = 0 --genome raw score/avarage number of genomes in species
  genome.network = {} --holds neurons
  genome.score = 0 --genome raw score
  genome.Fitness = 0
  genome.maxInnovation = 0
  genome.maxNeuron = 0
  genome.mutationChance = genomeMutationChance
  genome.flipSign = genomeFlipChance
  genome.activationChance = genomeActivationChance
  genome.weightDecrementChance = genomeDecrementChance
  genome.linkMutationChance = genomeLinkMutationChance
  genome.nodeMutationChance  = genomeNodeMutationChance
  genome.step = genomeStepValue
  return genome
end

 function newNeuron()
	local neuron = {}
  neuron.geneInnovationIndex = {} --stores innovation number of the connection gene its connected to
  neuron.inStatus = 0 --0 if an input 1 if an output (matches with node array index) 2 if normal neuron
  neuron.inputNumber = 0 --this is just for monitoring it on input nodes its not used anywhere else
	neuron.value = 0.0 -- current neuron value
	return neuron
end



 function BuildNetwork(genome)
  innovationNumber = 1
  maxNeuron = 1
  --add neurons
  for c=1, #TestOutputs do
    tempO = newNeuron()
    tempO.inStatus = 1 --signifies output neuron
    tempO.inputNumber = c
    genome.network[maxNeuron] = tempO
    maxNeuron = maxNeuron + 1
  end

	for i=1,#TestInputs do
    tempN = newNeuron()
    tempN.inStatus = 0 --signifies input neuron
    tempN.inputNumber = i
    for j = 1, #TestOutputs do
    tempConnectionGene = connectionGene()
    tempConnectionGene.input = maxNeuron
    tempConnectionGene.innovation = innovationNumber
    tempConnectionGene.out = j
    table.insert(genome.network[maxNeuron].geneInnovationIndex, tempConnectionGene.innovation)
    table.insert(genome.network[j].geneInnovationIndex, tempConnectionGene.innovation)
    if  genome.activationChance > math.random() then
      tempConnectionGene.status = true
     -- print("One or more of connection genes are enabled")
    else
      tempConnectionGene.status = false
     -- print("One or more of connection genes are disabled")
    end
    genome.genes[tempConnectionGene.innovation] = tempConnectionGene
    genome.maxInnovaion = tempConnectionGene.innovation
    innovationNumber = innovationNumber + 1
  end
  genome.network[maxNeuron] = tempN
  genome.maxNeuron = maxNeuron
  maxNeuron = maxNeuron + 1
	end
end

function selectionSort(genomes)
--sort by fitness
--int c = 0
for i = 1, #genomes do
maxFitnessIndex = i
for j = i + 1, #genomes do
--will now pick least fitness(closer to 100 weights sum) as highest fitness
if genomes[j].fitness < genomes[maxFitnessIndex].fitness then
maxFitnessIndex = j
end
end
tempg = genomes[i]
genomes[i] = genomes[maxFitnessIndex]
genomes[maxFitnessIndex] = tempg
end
end

function selectionSortForSpecies(speciesStore)
--speciesFitness
for i = 1, #speciesStore do
maxFitnessIndex = i
for j = i + 1, #speciesStore do
--(low to high)
--will now pick least fitness(closer to 100 weights sum) as highest fitness
if speciesStore[j].speciesFitness < speciesStore[maxFitnessIndex].speciesFitness then
maxFitnessIndex = j
end
end
tempg = speciesStore[i]
speciesStore[i] = speciesStore[maxFitnessIndex]
speciesStore[maxFitnessIndex] = tempg
end
end







function updateInputs(genome)
  inputCount = 1
  for i = 1, #genome.network do
    if genome.network[i].inStatus==0 and inputCount <= #TestInputs then
      --update input var
      genome.network[i].value = TestInputs[inputCount]
      inputCount = inputCount + 1
    end
  end
end





  function obtainOutputs(genome)
    fitness= 0
    sum = 0
    OutCount = 1
    value = 100
    --all weighs must add up to 100
    --print("first weight: "..genome.genes[1].weight)
  for i = 1, #genome.genes do
   --add weights
     sum = sum + genome.genes[i].weight
     --print("sum w loop:"..sum)
  end
  --print("net sum1: "..sum)
    --fitness = sum
    genome.score = (sum - 100)*-1
	if genome.score < 0 then
	genome.score = 100
	print("Too negative")
	end
	if genome.score > 100 then
	genome.score = 100
	print("Too positive")
	end
	--if #genome.network > MaxNodes then
	--genome.score = 100
	--end
    print("genome score :"..genome.score)
    --the lower the score, the better the genome
  --  print("net sum2: "..genome.score)
  end





function evaluateGenome(genome)
  updateInputs(genome)
  --print("NETWORK SIZE "..#genome.network)
  --print("NUMBER OF GENES IN NETWORK BEFORE EVALUATION "..#genome.genes)
  for i = 1, #genome.network do
    tempWout = {} --just the outs the ones we need(weights)
    oldValue = genome.network[i] --old neuron data
    for j = 1, #genome.network[i].weightIndex do
      tempWeightIndx = genome.network[i].weightIndex[j]
      for k = 1, #genome.genes do
        if tempWeightIndx == genome.genes[k].innovation and genome.genes[k].status == true and genome.genes[k].out == genome.network[i] then
          table.insert(tempWout,genome.genes[k])
        end
      end
    end
      sum = 0
      activation = 0
      for m = 1, #tempWout do
      sum = sum + (tempWout[m].input.value * tempWout[m].weight)
      end
      activation = sigmoid(sum)
      if #tempWout~=0 then
        genome.network[i].value = activation
      for n = 1, #tempWout do
        tempWout[n].out.value = activation
      end
        end
  end
  obtainOutputs(genome)
  --print("NUMBER OF GENES IN NETWORK AFTER EVALUATION "..#genome.genes)
end





function createStartingPopulation(number)
  genomesCreated = {}
  for i = 1, number do
    testPop = createNewGenome()
    BuildNetwork(testPop)
    if testPop.mutationChance > math.random() then
      mutateNodeGene(testPop)
    end
    if testPop.mutationChance > math.random() then
      mutateConnectionGene(testPop)
    end
    --evaluateNetwork(testPop)
    table.insert(genomesCreated,testPop)
   -- print("population no: "..i)
  end
  return genomesCreated
end




function createNewSpecies()
speciesMap = {}
speciesMap.genomes = {}
speciesMap.genomeMascot = createNewGenome()
--speciesMap.overallFitness = 0 --adjusted fitness(accumilation of fitness of all genes)
speciesMap.speciesFitness = 0
--speciesMap.attatchedSpecie = 0 --keeps track of where this species is in the species store
return speciesMap
end





function getAverageWeightDifference(genome,mascot)
  weightSum = 0
  for i = 1,#genome.genes do
    weightSum = weightSum + genome.genes[i].weight
  end
  for i = 1,#mascot.genes do
    weightSum = weightSum + mascot.genes[i].weight
  end
  Average = weightSum/(#genome.genes + #mascot.genes)
  return Average
end






function pickMascot(species)
  species.genomeMascot = species.genomes[math.random(1,#species.genomes)]
end



function resetSpecies(species)
  species.speciesFitness = 0
  species.genomeMascot = createNewGenome()
  species.genomes = {}
end



function calculateSingleSpeciesFitness(species) --evaluates a genome and gives genome fitness
for i=1,#species.genomes do
  --print("score in fit: "..species.genomes[i].score)
  species.genomes[i].fitness = species.genomes[i].score/#species.genomes
  species.speciesFitness = species.speciesFitness + species.genomes[i].fitness
  end
end




function calculateFitnessOfAllSpecies(species)
for i = 1, #speciesStore do
calculateSingleSpeciesFitness(speciesStore[i])
end
end



function bestGenomesInAllSpecies(species)

for i = 1, #species do
bestGenomesForNextGeneration(species[i])
end
end



--looks at best genomes in a single species and cuts out the rest
function bestGenomesForNextGeneration(specie)
  --pick best/fittest genomes in each species for next generation(top 3)
  culledGenomes = {}
  selectionSort(specie.genomes)
--top 3 fitness(lower fitness == better according to my current fitness func)
  for i = 1,#specie.genomes do
    if i == 10 then
    break
  end
  table.insert(culledGenomes,specie.genomes[i]) --reworks genomes in each specie
  table.insert(InitialPopulation1,specie.genomes[i]) --each genome selected is stored for next gen breeding
  end
specie.genomes = culledGenomes
print("First 10 genomes picked")
end

function selectBiasedSpeciesRandom()
--returns a species selected at random. Species with higher fitness is favoured(probability stuff)
chosenSpecie = {}
completeFitness = 0
for i = 1, #speciesStore do
completeFitness = completeFitness + speciesStore[i].speciesFitness
if next(speciesStore[i].genomes) == nil then
--print("This species has no genomes in select biased. This either means, 1) there was a species created but the genomes were not added")
end
end

print("CPS: "..completeFitness)

randomSelect = math.random()*completeFitness

print("CPSRSV: "..randomSelect)


countFit = 0

for i = 1, #speciesStore do
countFit = countFit + speciesStore[i].speciesFitness
--PROBABILITY SHIT (I FLIPPED THE SIGN. SPECIES WITH A LOWER FITNESS SCORE ARE MORE LIKEL TO BE SELECTED. THIS SIGN WILL BE FLIPPED DEPENDING ON FITNESS FUNC)
if countFit <= randomSelect then
chosenSpecie = speciesStore[i]
print("SPECIES CHOSEN")
break
end
end
if #speciesStore < 2 then
chosenSpecie = speciesStore[1]
end
if next(chosenSpecie) == nil then
print("Chosen species is empty. No species was able to meet ther random threshold")
end
return chosenSpecie
end

function selectRandomBiasedGenome(specie)
chosenGenome = {}
completeFitness = 0
for i = 1, #specie.genomes do
completeFitness = completeFitness + specie.genomes[i].fitness
end

print("CPG: "..completeFitness)

randomSelect = math.random()*completeFitness

print("CPGRSV: "..randomSelect)

countFit = 0

for i = 1, #specie.genomes do
countFit = countFit + specie.genomes[i].fitness
--PROBABILITY SHIT (I FLIPPED THE SIGN. GENOMES WITH A LOWER FITNESS SCORE ARE MORE LIKEL TO BE SELECTED. THIS SIGN( < OR >) WILL BE FLIPPED DEPENDING ON FITNESS FUNC)
if countFit <= randomSelect then
chosenGenome = specie.genomes[i]
print("GENOME CHOSEN")
print(#specie.genomes[i].genes)
break
end
end

return chosenGenome
end


--breed next generation from parent genomes
function nextGenerationMaker(Genomes)
 -- offSprings = {}
 print("number of genomes: "..#Genomes)
 --I GET IT NOW...IT DOESNT BREED OR ANYTHING BECAUSE NUMBER OF SPECIES > = POP SIZE
 --print("pop size: "..initialPopulationSize)
for i = #Genomes,initialPopulationSize do
  s1 = selectBiasedSpeciesRandom()
  if next(s1) == nil then
  print("RECOGNIZED")
  end
  if next(s1)~=nil then
  g1 = selectRandomBiasedGenome(s1)
  print("picked p1")
  else
  g1 = {}
  print("Species didnt meet threshold")
  end
  print("G1 DONE")
  if next(s1)~=nil then
  g2 = selectRandomBiasedGenome(s1)
  print("picked p2")
  else
  print("Species didnt meet threshold")
  end
--print("g1 network: "..#g1.network)
--print("g2 network: "..#g2.network)
if next(g1) == nil or next(g2) == nil then
print("One of us is barren because the expected random multiple threshold was not met :(")
else
  print("breeding............")
  if g1.fitness < g2.fitness then
  g3 = crossover(g1,g2)
  else
  g3 = crossover(g2,g1)
  end
  g3 = crossover(g1,g2)
  pointMutateGenome(g3)
  if g3.mutationChance > math.random() then
   if g3.nodeMutationChance > math.random() then
    mutateConnectionGene(g3)
  end
  if g3.linkMutationChance > math.random() then
    mutateNodeGene(g3)
  end
  end
  table.insert(Genomes,g3)
end
end
--print("pop size after breeding: "..#Genomes)
--return offSprings
end



function initializeSpecies(genome)
  newSpecies = createNewSpecies()
  table.insert(newSpecies.genomes,genome)
  newSpecies.mascot = genome
  table.insert(speciesStore,newSpecies)
end
--a random gene is chosen from a population
--each gene in the population is measured against the mascot
--set of genomes has to be from previous generation








function generateSpecies(Genomes)
  for i = 1,#Genomes do
    --loop through species
    --first check if species is empty or not
    if #speciesStore ~=0  then
      --print("species store contains speceis")
    for j = 1, #speciesStore do
      --choose mascot randomly
      pickMascot(speciesStore[j])
      s = speciationValue(Genomes[i],speciesStore[j].mascot)
      --if similar, add to mascot's species genomes
      if s < speciesDistance then
          table.insert(speciesStore[j].genomes,Genomes[i])
      break
      end
      --if they are too different, crate another species and put genome as mascot
      if s > speciesDistance then
      initializeSpecies(Genomes[i])
      break
      end
    end
  else
    --if species is nill, create one
    --print("no species in species store")
    initializeSpecies(Genomes[i])
    end
  end
  --print("Species foind: "..#speciesStore)
  end




function speciationValue(genome,mascot)
  x,y = disjointGenes(genome,mascot)
  e,f = excessGenes(genome,mascot)
  avg = getAverageWeightDifference(genome,mascot)
  a = #genome.genes
  b = #mascot.genes
  c= 0
  if a > b then
    c =a
  else
    c = b
    end
  speciationValueRes = (disjointmentConstant*#x/c)+(excessGenesConstant*#e/c)+(weightImportanceConstant*avg)
  --print("speciation value "..speciationValueRes)
  return speciationValueRes
  end
--general formula on paper: ( c1*E/N + c2*D/N + c3*W)



--takes entire population and calculates score of each genome
function evaluatePopulation(Population)
for i = 1, #Population do
evaluateGenome(Population[i])
end
end

--generateSpecies(Genomes)
--cullSpecies
--eliminates the last 1/4 of species and preserves the rest
function cullSpecies(speciesStore)
--1/4*(n + 1) = nth first quartile term
--so, I wanna preserve 3/4, delete the 1/4
selectionSortForSpecies(speciesStore) --sort species first before culling
culledSpecies = {}
Q3 = 3/4 * (#speciesStore + 1)
if Q3 < 1 then
Q3 = 1
end
for i = 1, # speciesStore do
if i <=Q3 then
table.insert(culledSpecies,speciesStore[i])
--print("Specie picked")
end

if next(speciesStore[i].genomes) == nil then
print("This species has no genomes. This either means, 1) there was a species created but the genomes were not added")
end


end
speciesStore = culledSpecies
end




InitialPopulation = createStartingPopulation(initialPopulationSize)

----------------------V2---------------------------------------------------------------------

for i = 1,NumberOfGenerations do
print("GENERATION: "..i)
evaluatePopulation(InitialPopulation)
generateSpecies(InitialPopulation)
calculateFitnessOfAllSpecies(speciesStore)
cullSpecies(speciesStore)
bestGenomesInAllSpecies(speciesStore) --culls out weak genomes in each species and selects top 3, best genomes stored for next gen
nextGenerationMaker(InitialPopulation1)
InitialPopulation = InitialPopulation1
InitialPopulation1 = {}
print("NUMBER OF SPECIES: "..#speciesStore)
if i~=NumberOfGenerations then
speciesStore = {}
end
end

print("NUMBER OF SPECEIES: "..#speciesStore)
print("SPECIES WITH HIGHEST FITNESS: "..speciesStore[1].speciesFitness)
for j = 1,#speciesStore do
for  k = 1, #speciesStore[j].genomes do
for l = 1, #speciesStore[j].genomes[k].genes do
print("WEIGHT VALUE: "..speciesStore[j].genomes[k].genes[l].weight)
end
end
end
