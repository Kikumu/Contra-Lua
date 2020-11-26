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


local TestInputs = {}
TestInputs[1] = 1
TestInputs[2] = 3
TestInputs[3] = 4
TestInputs[4] = 6
TestInputs[5] = 7

local TestOutputs = {}
TestOutputs[1] = 0.8
TestOutputs[2] = 0.8
TestOutputs[3] = 0.8
TestOutputs[4] = 0.8

local TestOutputs1 = {}
TestOutputs1[1] = 0.50
TestOutputs1[2] = 0.20
TestOutputs1[3] = 0.440
TestOutputs1[4] = 0.100

--stores all species
local speciesStore = {}

--used to create species
--GEN --SPECIES--GENOME


--geneActivationChance = 0.3

genomeFlipChance = 0.2
genomeDecrementChance = 0.2
genomeMutationChance = 0.3
genomeActivationChance = 0.25
genomeLinkMutationChance = 0.5


--genomePointMutateChance = 0.4
selectionChance = 0.3 --in use
initialPopulationSize = 100 --in use

crossOverChance = 0.75
NumberOfGenerations = 10000

--species classification variables

disjointmentConstant = 0.2 --c1
excessGenesConstant = 0.3 --c2
weightImportanceConstant = 0.1 --c3
speciesDistance = 0.06 --in use
MaxNodes = 1000

--E is number of disjointed connection genes
--D is number of excess connection genes
--W is weight value
--N number of connection genes
function softMax(Outputs)
  --val = val/sum of exp(val) entire val output net including val itself
  end
function createNewSpecies()
speciesMap = {}
speciesMap.genomes = {}
speciesMap.genomeMascot = createNewGenome()
--speciesMap.overallFitness = 0 --adjusted fitness(accumilation of fitness of all genes)
speciesMap.speciesFitness = 0
speciesMap.attatchedSpecie = 0 --keeps track of where this species is in the species store
return speciesMap
end

function pickMascot(species)
  species.genomeMascot = species.genomes[math.random(1,#species.genomes)]
  end

function resetSpecies(species)
  species.overallFitness = 0
  species.genomeMascot = species.species[math.random(1,#species.genomes)]
  species.genomes = {}
  end

function calculateFitness() --evaluates a genome and gives genome fitness
  --from best genomes chosen from next generation, if not enough to fit poopulation of next generation, thats when crossover and breeding comes through
end

function fitnessComparator(genome1,genome2)
  --if genome is greater print one else print/return 0
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

--a random gene is chosen from a population
--each gene in the population is measured against the mascot
--set of genomes has to be from previous generation
function generateSpecies(setOfGenomes)
  c = 0
  for i = 1,#setOfGenomes do
    --place genome in a species and move on to the next ay
    if #speciesStore > 0 then
    for j = 1, #speciesStore do
      pickMascot(speciesStore[j])
      print("species : "..j)
      print("species genome numbers: "..#speciesStore[j].genomes)
      print("number of species: "..#speciesStore)
      s = speciationValue(setOfGenomes[i],speciesStore[j].mascot)
      --print("speciation value "..s)
      if s < speciesDistance then --0.04
          --speciesStore[j].speciesFitness = (speciesStore[j].speciesFitness + setOfGenomes[i].fitness)/(#speciesStore[j].genomes + 1)
          speciesStore[j].speciesFitness = (speciesStore[j].speciesFitness + setOfGenomes[i].fitness)
          table.insert(speciesStore[j].genomes,setOfGenomes[i])
          --print("here1")
      break
    else
    newSpecies = createNewSpecies()
    table.insert(newSpecies.genomes,setOfGenomes[i])
    newSpecies.mascot = setOfGenomes[i]
    table.insert(speciesStore,newSpecies)
    --speciesStore[#speciesStore].speciesFitness = (speciesStore[#speciesStore].speciesFitness + setOfGenomes[i].fitness)/(speciesStore[#speciesStore].genomes)
    speciesStore[#speciesStore].speciesFitness = (speciesStore[#speciesStore].speciesFitness + setOfGenomes[i].fitness)
    --print("here2")
      break
      end
    end
  else
    --if theres nothing in the store, create a new species
    newSpecies = createNewSpecies()
    table.insert(newSpecies.genomes,setOfGenomes[i])
    newSpecies.mascot = setOfGenomes[i]
    table.insert(speciesStore,newSpecies)
    --speciesStore[#speciesStore].speciesFitness = (speciesStore[#speciesStore].speciesFitness + setOfGenomes[i].fitness)/(speciesStore[#speciesStore].genomes)
    speciesStore[#speciesStore].speciesFitness = (speciesStore[#speciesStore].speciesFitness + setOfGenomes[i].fitness)
    --print("here3")
    end
  end
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
  print("speciation value "..speciationValueRes)
  return speciationValueRes
  end
--general formula on paper: ( c1*E/N + c2*D/N + c3*W)


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
print("A GENE HAS BEEN MUTATED BY NODE")
end

function pointMutateGenome(genome)
  for i = 1, #genome.genes do
    --activation chance
    --weight inc/dec probability
    --flip sign probability
    if genome.flipSign > math.random() then
      genome.genes[i].weight = -genome.genes[i].weight
    end
    if genome.weightDecrementChance > math.random() then
     genome.genes[i].weight = genome.genes[i].weight - (genome.genes[i].weight * genome.step)
   else
      genome.genes[i].weight = genome.genes[i].weight * genome.step
    end

    if genome.activationChance > math.random() then
      genome.genes[i].status = true
    else
      genome.genes[i].status = false
    end
    end
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
function selectionSort(DisjointedGenesArr, ExcessGenesArr, MatchingGenesArr)
--genome = createNewGenome()
genesToSort = {}
accumilateGenesForSorting(DisjointedGenesArr,genesToSort)
accumilateGenesForSorting(ExcessGenesArr,genesToSort)
accumilateGenesForSorting(MatchingGenesArr,genesToSort)
--genome.genes = genesToSort
return genesToSort
end



function BuildNetworkOfChildGene(gene)
  --build i and o
  neurons2 = {}
  genome = createNewGenome()
  --connect/create inputs
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
        --print("input detected")
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
         --print("output detected which shouldnt be done(illegal)")
         gene[j].out = TempOutput
        table.insert(TempOutput.weightIndex,gene[j].innovation)
      end
      if gene[j].out.inStatus == 1 and gene[j].out.inputNumber==o then
        TempOutput.inStatus = 1
        TempOutput.inputNumber = o
         --print("output detected is allowed")
         gene[j].out = TempOutput
        table.insert(TempOutput.weightIndex,gene[j].innovation)
      end
    end
    table.insert(genome.network,TempOutput)
  end
  print("number of neurons initially"..#neurons2)
  for k = 1,#gene do
    if gene[k].input.inStatus == 2 then
        table.insert(neurons2,gene[k].input)
        --make gene notice this new point in neurons2
        gene[k].input = neurons2[#neurons2]
       print("normal neuron detected 1")
      end
      --if out connection gene is a normal neuron
      if gene[k].out.inStatus == 2 then
        table.insert(neurons2,gene[k].out)
        gene[k].out = neurons2[#neurons2]
        print("normal neuron detected 2")
      end
    end
  print("number of neurons after rolling in genes"..#neurons2)
  --create and connect other neurons
  --other neurons are labelled "2"
  --they can either be "in" connection gene or "out"
  --how many of "2"'s are there? --#neurons2
  --how many are duplicates of the other?
  isDuplicate = 0
  --indexArr = {}
  for d = 1, #neurons2 do
    isDuplicate = 0
    for d1 = d+1, #neurons2 do
      if neurons2[d] == neurons2[d1] then
       -- isDuplicate = 1
        table.remove(neurons2,d1)
        --table.insert(indexArr,d1)
      end
    end

  end
  print("number of neurons after clean  "..#neurons2)
    for i = 1, #neurons2 do
      table.insert(genome.network,neurons2[i])
    end
    genome.genes = gene
    --genome.genes = copyGene(gene)
    return genome
  end
--function newNeuron()
--	local neuron = {}
--  neuron.weightIndex = {} --stores innovation number of the connection gene its connected to
--  neuron.inStatus = 2 --0 if an input 1 if an output (matches with node array index) 2 if normal neuron
--  neuron.inputNumber = 0 --this is just for monitoring it on input nodes its not used anywhere else
--  neuron.outputNumber = 0 --this is just for monitoring it on output nodes its not used anywhere else
--	neuron.incoming = {} --data from previous thingis
--	neuron.value = 0.0 -- current neuron value, is alwas a zero
--	return neuron
--end
--breed
function crossover(genome1,genome2) --genome in the sense that you are passing to this function a set of genomes
m = matchingGenes(genome1,genome2)
--print("number of matching neurons: "..#n1)
print("number of matching genes: "..#m)
d = disjointGenes(genome1,genome2)
--print("number of disjoint neurons: "..#n2)
print("number of disjoint genes: "..#d)
e = excessGenes(genome1,genome2)
--print("number of excess neurons: "..#n3)
print("number of excess genes: "..#e)
--extract neurons from genes and put in a network
genes = selectionSort(d, e, m)
--genome.network = nT1
print("number of all genes: "..#genes)
return BuildNetworkOfChildGene(genes)
end

function species(genes)
speciesPool = {}
return speciesPool
end

--neurons should match with innovation number
function createNewGenome()
  genome = {} --holds the gene information network (general info)
  genome.genes = {} --weight info (connection genes info)
  genome.fitness = 0
  genome.network = {} --holds neurons
  genome.score = 0
  genome.mutationChance = math.random() --genomeMutationChance
  genome.flipSign = math.random() --genomeFlipChance
  genome.activationChance = math.random() -- genomeActivationChance
  genome.weightDecrementChance = math.random() --genomeDecrementChance
  genome.linkMutationChance = math.random() --genomeLinkMutationChance
  genome.step = math.random(0.5,1.5)
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

--for testing purposes
  function obtainOutputs(genome)
    OutCount = 1
    fitness= 0
  for i = 1, #genome.network do
    if genome.network[i].inStatus==1 and OutCount <= #TestOutputs then
      --update input var
      --fitness = genome.network[i].value - TestOutputs[OutCount]
      print("New Output"..genome.network[i].value)
      OutCount = OutCount + 1
      end
    end
    end

function evaluateGenome(genome)
  --obtain/update input neurons (input neurons have a state of 1)
  updateInputs(genome)
  --obtain all connections to a node and shit out output
  for i = 1, #genome.network do
    --print("old neuron value".. genome.network[i].value)
    tempWout = {} --just the outs the ones we need(weights)
    oldValue = genome.network[i] --old neuron data
    --update node AND gene
    --obtain node value
    --TempNode = genome.network[i]
    --print("genome network neuron old::  " ..genome.network[i].value)
    --print("weight index values " ..#genome.network[i].weightIndex) --how many innovation numbers are in here
    --obtain weights through index
    for j = 1, #genome.network[i].weightIndex do --loops through stored gene innovation numbers of neuron
      tempWeightIndx = genome.network[i].weightIndex[j] --current innovation number stored in neuron
      --print("innovation number to look for"..tempWeightIndx)
      for k = 1, #genome.genes do --look for neuron inn number in genome genes
        --if innovation number of neuron is equal to innovation number of connection gene in genome
        --and if this is not an input(remember, input is 0 normal is 2 output is 1)
        if tempWeightIndx == genome.genes[k].innovation and genome.genes[k].status == true and genome.genes[k].out == genome.network[i] then
          --add to tempW
         --print("Associated gene in"..genome.genes[k].input.value)
         --print("Associated WEIGHT "..genome.genes[k].weight)
         --print("Associated gene out"..genome.genes[k].out.value)
         --add to tempW(all ins and outs connection genes)
          table.insert(tempWout,genome.genes[k])--this will now hold all associated connection genes
        end
      end
    end
      --after finding all associated genes, filter ins and outs(if gene out == genome.network) and store temporarily
     -- print("found refined connection genes : "..#tempWout)
    --loop through all cases of tempWout, obtain ins, multiply by weights get new out value and REPLACE the gene with the new tempWouts
      sum = 0
      activation = 0
      --calculate sum
      for m = 1, #tempWout do
    -- print("input value from filtered "..tempWout[m].input.value)
    -- print("output value from filtered "..tempWout[m].out.value)
    -- print("weight value from filtered "..tempWout[m].weight)
      sum = sum + (tempWout[m].input.value * tempWout[m].weight)
     --print("sum: "..sum)
      end
      activation = sigmoid(sum)
      --print("activated")
    --print("activation value: "..activation)
      if #tempWout~=0 then
        genome.network[i].value = activation
      --replace all values of #tempwout.out with new activation value
      for n = 1, #tempWout do
        tempWout[n].out.value = activation
      end
        end

      --replace in actual gene as well (both ins and outs)
      --REPLACE NODE
    --print("new neuron value".. genome.network[i].value)
    --print("new genome i/o stats".. genome.network[i].inStatus)
  end
  pointMutateGenome(genome)
  obtainOutputs(genome)
  if genome.mutationChance > math.random() then
    mutateConnectionGene(genome)
  end
  if genome.linkMutationChance > math.random() then
    mutateNodeGene(genome)
  end

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

InitialPopulation = createStartingPopulation(initialPopulationSize)
--evaluate each genome in population and group into species
generateSpecies(InitialPopulation)
for i = 1,#speciesStore do
  print("Number of genomes in this species: "..#speciesStore[i].genomes)
  end
