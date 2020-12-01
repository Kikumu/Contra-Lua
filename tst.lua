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
local InitialPopulation = {}
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
print("Species first :"..#speciesStore)
--used to create species
--GEN --SPECIES--GENOME
--geneActivationChance = 0.3
genomeFlipChance = 0.20
genomeDecrementChance = 0.35
genomeMutationChance = 0.5
genomeActivationChance = 0.25
genomeLinkMutationChance = 0.5
genomeNodeMutationChance = 0.35
genomeStepValue = 0.85

--genomePointMutateChance = 0.4
selectionChance = 0.03 --in use
initialPopulationSize = 150 --in use
crossOverChance = 0.75
NumberOfGenerations = 10
--species classification variables
disjointmentConstant = 0.25 --c1
excessGenesConstant = 0.3 --c2
weightImportanceConstant = 0.35 --c3
speciesDistance = 0.30 --in use
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



--------------------------------------------------------------------

-----------------------------------------------------------------







--activation function
function sigmoid(x)
    return 1.0 / (1.0 + math.exp(-x))
end






function connectionGene()
local gene = {}
gene.input = 0
gene.out = 0
gene.weight = math.random()
gene.status = true
gene.innovation = 0 --ancestry monitor
return gene
end







--tested(needs reviewing) --TESTED AND REWORKED
function mutateConnectionGene(genome)

max = retunMaxInnovation(genome)
--grab 2 random nodes from network
v1 = math.random(1,#genome.network)
v2 = math.random(1,#genome.network)
--print("max val in connection gene"..max)
node1 = genome.network[v1] --1
node2 = genome.network[v2] --1
--loop through the 4 neurons and check if same innov number
state1 = 0
state2 = 0
state3 = 0
 --check if same innovation indexes
for i = 1, #node1.weightIndex do
  for j = 1, #node2.weightIndex do
    if node1.weightIndex[i] == node2.weightIndex[j] then
      state1 = 1
      end
  end
end
if state1 == 0 then
  if node1.inStatus~= 1 then --if input is not a network output node
    state2 = 1
    end
end
if state2 == 1 then
  if node2.inStatus~= 0 then --if input is not a network input node
    state3 = 1
    end
end
if state3 == 1 then
connectionGeneMutate1 = connectionGene()
connectionGeneMutate1.input = node1
connectionGeneMutate1.innovation = max + 1
--print("mutate a: ",connectionGeneMutate1.input.value)
table.insert(connectionGeneMutate1.input.weightIndex,connectionGeneMutate1.innovation) --register additional gene to neuron
genome.network[v1] = connectionGeneMutate1.input
connectionGeneMutate1.out = node2
table.insert(connectionGeneMutate1.out.weightIndex,connectionGeneMutate1.innovation) --register additional gene to neuron
genome.network[v2] = connectionGeneMutate1.out
--print("mutate b: ",connectionGeneMutate1.out.value)
print("a gene has been mutated by connection")
table.insert(genome.genes,connectionGeneMutate1)
end

--return genome
end





function retunMaxInnovation(genome)
maxInnovation = genome.genes[1].innovation
for i = 1, #genome.genes do
  if genome.genes[i].innovation > maxInnovation then
  maxInnovation = genome.genes[i].innovation
  end
end
return maxInnovation
end






--FULLY REWORKED
function mutateNodeGene(genome)
tstval = #genome.genes
maxInnovation = retunMaxInnovation(genome)
--print("Mutation max prev innovation value: "..maxInnovation)
genepos = math.random(1,#genome.genes)
connectionGeneTemp = genome.genes[genepos]
connectionGeneTemp.status = false
genome.genes[genepos] = connectionGeneTemp
connectionGeneMutate1 = connectionGene()
connectionGeneMutate1.input = connectionGeneTemp.input --attatch neuron
connectionGeneMutate1.innovation = maxInnovation + 1
table.insert(connectionGeneMutate1.input.weightIndex,connectionGeneMutate1.innovation)
tempNeuron = newNeuron()
tempNeuron.value = math.random()
table.insert(tempNeuron.weightIndex,connectionGeneMutate1.innovation)
connectionGeneMutate1.out = tempNeuron
table.insert(genome.genes,connectionGeneMutate1)
connectionGeneMutate2 = connectionGene()
connectionGeneMutate2.input = tempNeuron
connectionGeneMutate2.innovation = connectionGeneMutate1.innovation + 1
table.insert(tempNeuron.weightIndex,connectionGeneMutate2.innovation)
connectionGeneMutate2.out = connectionGeneTemp.out
table.insert(connectionGeneMutate2.out.weightIndex,connectionGeneMutate2.innovation)
table.insert(genome.network,tempNeuron)
table.insert(genome.genes,connectionGeneMutate2)
--print("A GENE HAS BEEN MUTATED BY NODE")
end






function pointMutateGenome(genome)
--mutateState = 0
print("S "..#genome.genes)
print("F"..genome.flipSign)
print("W"..genome.weightDecrementChance)
print("A"..genome.activationChance)

print("number of genes in pointMutate Before: "..#genome.genes)
  for i = 1, #genome.genes do
    if genome.flipSign > math.random() then
      genome.genes[i].weight = -genome.genes[i].weight
	  print("Sign flipped")
      break
    end
    if genome.weightDecrementChance > math.random() then
	 genome.genes[i].weight = genome.genes[i].weight - (genome.genes[i].weight * genome.step)
	  print("Weight decrement")
      break

   else
	   genome.genes[i].weight = genome.genes[i].weight + (genome.genes[i].weight * genome.step)
	 print("Weight increment")
     break
    end

    if genome.activationChance > math.random() then
      genome.genes[i].status = true
	  print("gene enabled")
      break
    else
      genome.genes[i].status = false
	   print("gene disabled")
      break
    end
  end
print("number of genes in pointMutate After: "..#genome.genes)
end




--for crossover purposes(Goal is to add to g1)tested--TESTED AND REWORKED
function matchingGenes(genome1,genome2)
  --genes in which innovations match
  matchedGenes = {}
  for i = 1, #genome1.genes do
    for j = 1, #genome2.genes do
      if genome1.genes[i].innovation == genome2.genes[j].innovation then
        if selectionChance > math.random() then
           table.insert(matchedGenes,genome1.genes[i])
        else
          table.insert(matchedGenes,genome2.genes[j])
        end
      end
    end
  end
  return matchedGenes
end
--(takes disjointed genes from g2), goal is to add to G1
--the genes in the middle





function disjointGenes(genome1,genome2)
disjointedGenes = {}
max = retunMaxInnovation(genome1)
max2 = retunMaxInnovation(genome2)
--print("max in genome1 "..max)
--print("max in genome2 "..max2)
found = 0
excess = 0
for i = 1, #genome1.genes do
  found = 0
  excess = 0
  for j = 1, #genome2.genes do
    if genome2.genes[j].innovation <= max and genome1.genes[i].innovation == genome2.genes[j].innovation then
      found = 1
    end
    if genome2.genes[j].innovation  > max and found == 0 then
      excess = 1
    end
  end
  if found == 0 and excess == 0 then
    table.insert(disjointedGenes,genome2.genes[i])
  end
end
return disjointedGenes
end





function excessGenes(genome1,genome2)
ExcessGenesTable = {}
max = retunMaxInnovation(genome1)

  for j = 1, #genome2.genes do
    if genome2.genes[j].innovation > max then
      table.insert(ExcessGenesTable,genome2.genes[j])
      end
  end
return ExcessGenesTable
end






function accumilateGenesForSorting(genes,tbs)
  for i = 1, #genes do
    table.insert(tbs,genes[i])
  end
  return tbs
end






--sorted child gene according to innovation number
function arrayCombinerForCrossover(DisjointedGenesArr, ExcessGenesArr, MatchingGenesArr)
--genome = createNewGenome()
genesToSort = {}
accumilateGenesForSorting(DisjointedGenesArr,genesToSort)
accumilateGenesForSorting(ExcessGenesArr,genesToSort)
accumilateGenesForSorting(MatchingGenesArr,genesToSort)
--genome.genes = genesToSort
return genesToSort
end









function crossover(genome1,genome2) --genome in the sense that you are passing to this function a set of genomes
matched = matchingGenes(genome1,genome2)
--print("number of matching neurons: "..#n1)
print("number of matching genes: "..#matched)
disjointed = disjointGenes(genome1,genome2)
--print("number of disjoint neurons: "..#n2)
print("number of disjoint genes: "..#disjointed)
excess = excessGenes(genome1,genome2)
--print("number of excess neurons: "..#n3)
print("number of excess genes: "..#excess)
--extract neurons from genes and put in a network
genes = arrayCombinerForCrossover(matched,disjointed,excess)
--genome.network = nT1
--print("number of all genes: "..#genes)

genome = BuildNetworkOfChildGene(genes)
print("Crossover genes number before: "..#genome.genes)

if genome.mutationChance > math.random() then
print("Point mutate attempt")
pointMutateGenome(genome)
--print("pintMutate")
end
if genome.linkMutationChance > math.random() then
print("ConnectMutate Attempt")
mutateConnectionGene(genome)
end
if genome.nodeMutationChance > math.random() then
print("NodeMutate Attempt")
mutateNodeGene(genome)
end

print("Crossover genes number after mutations: "..#genome.genes)
return genome
end






function BuildNetworkOfChildGene(gene)
  --build i and o
  neurons2 = {}
  genome = createNewGenome()
  --connect/create inputs
  print ("number of genes in child before in BuildChildNetwork "..#gene)
  for i = 1, #TestInputs do
    TempInput = newNeuron()
    --how many in gene are connected to this one?
    --how to check,check if gene is an input (a zero by calling instats) and by checking the input number
    --if they are equal, add the connection gene to this neuron
    for j = 1, #gene do
      if gene[j].input.inStatus == 0 and gene[j].input.inputNumber==i then
        TempInput.inStatus = gene[j].input.inStatus
        TempInput.inputNumber = i
        table.insert(TempInput.weightIndex,gene[j].innovation)
        gene[j].input = TempInput
      --  print("input detected")
      end
    end
    table.insert(genome.network,TempInput)
end
--
--connect outputs
  for o = 1, #TestOutputs do
    TempOutput = newNeuron()
    for j = 1, #gene do
      if gene[j].input.inStatus == 1 and gene[j].input.inputNumber==o then
        TempOutput.inStatus = 1
        TempOutput.inputNumber = o
       --  print("output detected which shouldnt be done(illegal)")
         gene[j].out = TempOutput
        table.insert(TempOutput.weightIndex,gene[j].innovation)
      end
      if gene[j].out.inStatus == 1 and gene[j].out.inputNumber==o then
        TempOutput.inStatus = 1
        TempOutput.inputNumber = o
      --   print("output detected is allowed")
         gene[j].out = TempOutput
        table.insert(TempOutput.weightIndex,gene[j].innovation)
      end
    end
    table.insert(genome.network,TempOutput)
  end
 -- print("number of neurons initially"..#neurons2)
 -- print("number of neurons in networ initially"..#genome.network)
  for k = 1,#gene do
    if gene[k].input.inStatus == 2 then
        table.insert(neurons2,gene[k].input)
        --make gene notice this new point in neurons2
        gene[k].input = neurons2[#neurons2]
    --   print("normal neuron detected 1")
      end
      --if out connection gene is a normal neuron
      if gene[k].out.inStatus == 2 then
        table.insert(neurons2,gene[k].out)
        gene[k].out = neurons2[#neurons2]
      --  print("normal neuron detected 2")
      end
    end


  for d = 1, #neurons2 do
    isDuplicate = 0
    for d1 = d+1, #neurons2 do
      if neurons2[d] == neurons2[d1] then
        table.remove(neurons2,d1)
      end
    end

  end
--  print("number of neurons after clean  "..#neurons2)
    for i = 1, #neurons2 do
      table.insert(genome.network,neurons2[i])
    end
    genome.genes = gene
    --genome.genes = copyGene(gene)
    print("number of genes in child after in BuildChildNetwork"..#genome.genes)
    return genome
  end








--neurons should match with innovation number
function createNewGenome()
  genome = {} --holds the gene information network (general info)
  genome.genes = {} --weight info (connection genes info)
  genome.fitness = 0 --genome raw score/avarage number of genomes in species
  genome.network = {} --holds neurons
  genome.score = 0 --genome raw score
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
  neuron.weightIndex = {} --stores innovation number of the connection gene its connected to
  neuron.inStatus = 2 --0 if an input 1 if an output (matches with node array index) 2 if normal neuron
  neuron.inputNumber = 0 --this is just for monitoring it on input nodes its not used anywhere else
	neuron.value = 0.0 -- current neuron value
	return neuron
end

--gene.input = 0
--gene.out = 0
--gene.weight = math.random()
--gene.status = true
--gene.innovation = 0 --ancestry monitor

--you only need to do this ONCE per genome initialization










 function BuildNetwork(genome)
  innovationNumber = 1
  --add neurons
  for c=1, #TestOutputs do
    tempO = newNeuron()
    tempO.value = TestOutputs[c]
    tempO.inStatus = 1 --signifies output neuron
    tempO.inputNumber = c
    table.insert(genome.network,tempO)
  end

	for i=1,#TestInputs do
    tempN = newNeuron()
    tempN.value = TestInputs[i]
    tempN.inStatus = 0 --signifies input neuron
    tempN.inputNumber = i
    --where you make a connection gene is where u need to put innovation numbers
    for j = 1, #TestOutputs do
    --create gene and link innovation number
    tempConnectionGene = connectionGene()
    tempConnectionGene.input = tempN
    tempConnectionGene.innovation = innovationNumber
    randN = math.random()
    --print("random value "..randN)
    if  genome.activationChance > math.random() then
      tempConnectionGene.status = true
     -- print("One or more of connection genes are enabled")
    else
      tempConnectionGene.status = false
     -- print("One or more of connection genes are disabled")
    end
    --call network
    tempConnectionGene.out = genome.network[j]
    --add connection gene
    table.insert(tempConnectionGene.out.weightIndex,innovationNumber)
    --reassign
    genome.network[j] = tempConnectionGene.out
    table.insert(tempN.weightIndex,innovationNumber)
    table.insert(genome.genes,tempConnectionGene)
    innovationNumber = innovationNumber + 1
    end
    table.insert(genome.network,tempN)
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
	if #genome.network > MaxNodes then
	genome.score = 100
	end
    print("genome score :"..genome.score)
    --the lower the score, the better the genome
  --  print("net sum2: "..genome.score)
  end





function evaluateGenome(genome)
  updateInputs(genome)
  print("NETWORK SIZE "..#genome.network)
  print("NUMBER OF GENES IN NETWORK BEFORE EVALUATION "..#genome.genes)
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
  print("NUMBER OF GENES IN NETWORK AFTER EVALUATION "..#genome.genes)
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
    if i == 3 then
    break
  end
  table.insert(culledGenomes,specie.genomes[i]) --reworks genomes in each specie
  table.insert(InitialPopulation1,specie.genomes[i]) --each genome selected is stored for next gen breeding
  end
specie.genomes = culledGenomes

end

function selectBiasedSpeciesRandom()
--returns a species selected at random. Species with higher fitness is favoured(probability stuff)
chosenSpecie = {}
completeFitness = 0
for i = 1, #speciesStore do
completeFitness = completeFitness + speciesStore[i].speciesFitness
end

randomSelect = math.random()*completeFitness

countFit = 0

for i = 1, #speciesStore do
countFit = countFit + speciesStore[i].speciesFitness
if countFit >= randomSelect then
chosenSpecie = speciesStore[i]
break
end
end
return chosenSpecie
end

function selectRandomBiasedGenome(specie)
chosenGenome = {}
completeFitness = 0
for i = 1, #specie.genomes do
completeFitness = completeFitness + specie.genomes[i].fitness
end

randomSelect = math.random()*completeFitness

countFit = 0

for i = 1, #specie.genomes do
countFit = countFit + specie.genomes[i].fitness
if countFit >= randomSelect then
chosenGenome = specie.genomes[i]
print("GENOME CHOSEN"..#chosenGenome.genes)
break
end
end
return chosenGenome
end


--breed next generation from parent genomes
function nextGenerationMaker(Genomes)
 -- offSprings = {}
 print("number of genomes: "..#Genomes)
 print("pop size: "..initialPopulationSize)
for i = #Genomes,initialPopulationSize do
  s1 = selectBiasedSpeciesRandom()
  g1 = selectRandomBiasedGenome(s1)
  g2 = selectRandomBiasedGenome(s1)
--  print("g1 network: "..#g1.network)
 -- print("g2 network: "..#g2.network)
if #g1.genes == 0 or #g2.genes == 0 then
print("One of us is barren for some reason :(")
else
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
print("pop size after breeding: "..#Genomes)
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
    print("no species in species store")
    initializeSpecies(Genomes[i])
    end
  end
  print("Species foind: "..#speciesStore)
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
end
end
speciesStore = culledSpecies
end




InitialPopulation = createStartingPopulation(initialPopulationSize)
--evaluate each genome in population and group into species
--generateSpecies(InitialPopulation)

--[[for i = 1, NumberOfGenerations do
evaluatePopulation(InitialPopulation)
generateSpecies(InitialPopulation) --auto saves to species store
calculateFitnessOfAllSpecies(speciesStore)
bestGenomesInAllSpecies(speciesStore) --saves in separate temp Pop
nextGenerationMaker(InitialPopulation1)
InitialPopulation = InitialPopulation1
InitialPopulation1 = {}
speciesStore = {}
end]]

----------------------V2

for i = 1,NumberOfGenerations do
evaluatePopulation(InitialPopulation)
generateSpecies(InitialPopulation)
calculateFitnessOfAllSpecies(speciesStore)
cullSpecies(speciesStore)
bestGenomesInAllSpecies(speciesStore) --culls out weak genomes in each species and selects top 3, best genomes stored for next gen
nextGenerationMaker(InitialPopulation1)
InitialPopulation = InitialPopulation1
InitialPopulation1 = {}
if i~=NumberOfGenerations then
speciesStore = {}
end
end

print("NUMBER OF SPECEIES: "..#speciesStore)
for j = 1,#speciesStore do
for  k = 1, #speciesStore[j].genomes do
for l = 1, #speciesStore[j].genomes[k].genes do
print("WEIGHT VALUE: "..speciesStore[j].genomes[k].genes[l].weight)
end
end
end
