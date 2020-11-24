 local TestInputs = {}
TestInputs[1] = 12
TestInputs[2] = 39
TestInputs[3] = 34
TestInputs[4] = 629
TestInputs[5] = 639

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

nodeMutationChance = 0.5  --in use
geneMutationChance = 0.2  --in use
geneActivationChance = 0.5
selectionChance = 0.3 --in use
initialPopulationSize = 1000 --in use
speciesDistance = 0.23 --in use
crossOverChance = 0.75
NumberOfGenerations = 10000

--species classification variables

disjointmentConstant = 0.2 --c1
excessGenesConstant = 0.3 --c2
weightImportanceConstant = 0.1 --c3
MaxNodes = 11

--E is number of disjointed connection genes
--D is number of excess connection genes
--W is weight value
--N number of connection genes

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

function speciation(genome,mascot)
  newSpecies1 = {}
  newSpecies2 = {}
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
  speciationValue = (disjointmentConstant*#x/c)+(excessGenesConstant*#e/c)+(weightImportanceConstant*avg)
  print("speciation value "..speciationValue)
  if speciationValue > speciesDistance then
    --table.insert(newSpecies1,genome)
    print("specie1")
  else
    --table.insert(newSpecies2,genome)
    print("specie2")
    end
  end
--general formula on paper: ( c1*E/N + c2*D/N + c3*W)
--[[
TO DO:
1) Take it out for a test spin to test out functions created --DONE
2) Test propagation --Done
3) Neurons holder for each network created --DONE
4) Find a way to keep track of the fitness of each genome --done
5) Create function to organise a group of genes into different speciesPool
6) Rework how the network breeds --done
7) Fix how genes mutate --done
8) Integrate to game
]]

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
print("max val in connection gene"..max)
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
print("mutate a: ",connectionGeneMutate1.input.value)
table.insert(connectionGeneMutate1.input.weightIndex,connectionGeneMutate1.innovation) --register additional gene to neuron
genome.network[v1] = connectionGeneMutate1.input
connectionGeneMutate1.out = node2
table.insert(connectionGeneMutate1.out.weightIndex,connectionGeneMutate1.innovation) --register additional gene to neuron
genome.network[v2] = connectionGeneMutate1.out
print("mutate b: ",connectionGeneMutate1.out.value)
print("a gene has been mutated by connection")
table.insert(genome.genes,connectionGeneMutate1)
end
--return genome
end

function retunMaxInnovation(gene)
maxInnovation = gene.genes[1].innovation
for i = 1, #gene.genes do
  if gene.genes[i].innovation > maxInnovation then
  maxInnovation = gene.genes[i].innovation
  end
end
return maxInnovation
end


--FULLY REWORKED
function mutateNodeGene(genome)
tstval = #genome.genes
maxInnovation = retunMaxInnovation(genome)
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
print("This gene has been mutated by node")
end

--for crossover purposes(Goal is to add to g1)tested--TESTED AND REWORKED
function matchingGenes(genome1,genome2)
matchedGenes ={}
neuronsExtracted = {}
g1Length = #genome1.genes
g2Length = #genome2.genes
for i = 1, g1Length do
for j = 1, g2Length do
if genome1.genes[i].innovation == genome2.genes[j].innovation then
if selectionChance > math.random() then
table.insert(matchedGenes,genome1.genes[i])
table.insert(neuronsExtracted,genome1.genes[i].input)
table.insert(neuronsExtracted,genome1.genes[i].out)
print("picked a "..j)
else
table.insert(matchedGenes,genome2.genes[j])
table.insert(neuronsExtracted,genome2.genes[j].input) --grab neuron in
table.insert(neuronsExtracted,genome2.genes[j].out) --grab neuron out
print("picked b "..j)
end
end
end
end
return matchedGenes, neuronsExtracted
end
--(takes disjointed genes from g2), goal is to add to G1
--the genes in the middle
function disjointGenes(genome1,genome2)
maxInnovation = 0
disjointedGenes = {}
neuronsExtracted = {}
g2MaxInnovation = retunMaxInnovation(genome2)
g1MaxInnovation = retunMaxInnovation(genome1)
g2Length = #genome2
g1Length = #genome1
found = 0 --checks if connection gene is found in the thing or not
addLimiter = 0 --prevents adding of excess genes 1 add 0 dont add
--finding disjointed genes in g2
for i = 1, #genome2 do
  found = 0
  addLimiter = 0
  for j = 1, #genome1 do
    if genome1.genes[j].innovation < g2MaxInnovation then
      --if they are equal, matching gene found which is not what we want so nothing is saved
      addLimiter = 1 --add to disjointed at end if nothing is found
    if genome2.genes[i].innovation == genome1.genes[j].innovation then
      found = 1
    end
  end
    end
   if found == 0 and  addLimiter == 1 then
    table.insert(disjointedGenes,genome2.genes[i])
    table.insert(neuronsExtracted, genome2.genes[i].input)
    table.insert(neuronsExtracted, genome2.genes[i].out)
  end
end
return disjointedGenes,neuronsExtracted
end

function excessGenes(genome1,genome2)
ExcessGenesTable = {}
neuronsExtracted = {}
--find excess genes in g2(which are located in g1)
maxInnovationg1 =retunMaxInnovation(genome1)
maxInnovationg2 = retunMaxInnovation(genome2)
g1Length = #genome1.genes
g2Length = #genome2.genes
if maxInnovationg1 > maxInnovationg2 then
for i = 1, g1Length do
if genome1.genes[i].innovation > maxInnovationg2 then
--assumes g1 is fitter gene
table.insert(ExcessGenesTable, genome1.genes[i])
table.insert(neuronsExtracted,genome1.genes[i].out)
table.insert(neuronsExtracted,genome1.genes[i].input)
end
end
return ExcessGenesTable, neuronsExtracted
end
if maxInnovationg2 > maxInnovationg1 then
for i = 1, g2Length do
if genome2.genes[i].innovation > maxInnovationg1 then
  --assumes g2 is fitter gene
table.insert(ExcessGenesTable, genome2.genes[i])
table.insert(neuronsExtracted,genome2.genes[i].out)
table.insert(neuronsExtracted,genome2.genes[i].input)
end
end
end
return ExcessGenesTable, neuronsExtracted
end

function accumilateGenesForSorting(genes,tbs)
  for i = 1, #genes do
    table.insert(tbs,genes[i])
  end
  return tbs
end
--sorted child gene according to innovation number
function selectionSort(DisjointedGenesArr, ExcessGenesArr, MatchingGenesArr)
genome = createNewGenome()
genesToSort = {}
accumilateGenesForSorting(DisjointedGenesArr,genesToSort)
accumilateGenesForSorting(ExcessGenesArr,genesToSort)
accumilateGenesForSorting(MatchingGenesArr,genesToSort)
for i = 1, #genesToSort do
minInnovationIndex = i
for j = i + 1, #genesToSort do
if genesToSort[j].innovation < genesToSort[minInnovationIndex].innovation then
minInnovationIndex = j
end
end
tempg = genesToSort[i]
genesToSort[i] = genesToSort[minInnovationIndex]
genesToSort[minInnovationIndex] = tempg
end
genome.genes = genesToSort
return genome
end

function combineNeurons(n1,n2,n3)
  network = {}
  for i = 1,#n1 do
    table.insert(network,n1[i])
  end
  for i = 1,#n2 do
    table.insert(network,n2[i])
  end
  for i = 1,#n3 do
    table.insert(network,n3[i])
  end
  return network
end

function removeDuplicateNeurons(network)
  cleanNeurons = {}
  for i=1, #network do
    c = 1
    for j = 1, #network do
      if network[i] == network[j] then
        c = c + 1
        end
      end
      if c < 3 then
          table.insert(cleanNeurons,network[i])
        end
    end
    return cleanNeurons
  end
--breed
function crossover(genome1,genome2) --genome in the sense that you are passing to this function a set of genomes
m,n1 = matchingGenes(genome1,genome2)
d,n2 = disjointGenes(genome1,genome2)
e,n3 = excessGenes(genome1,genome2)
genome = selectionSort(d, e, m)
nT = combineNeurons(n1,n2,n3)
nT1 = removeDuplicateNeurons(nT)
genome.network = nT1
return genome
end

function species(genes)
speciesPool = {}
return speciesPool
end

--neurons should match with innovation number
function createNewGenome()
  genomeCluster = {} --holds the gene information network (general info)
  genomeCluster.genes = {} --weight info (connection genes info)
  genomeCluster.fitness = 0
  genomeCluster.network = {} --holds neurons
  return genomeCluster
end

 function newNeuron()
	local neuron = {}
  neuron.weightIndex = {} --stores innovation number of the connection gene its connected to
  neuron.inStatus = 2 --0 if an input 1 if an output (matches with node array index) 2 if normal neuron
  neuron.inputNumber = 0 --this is just for monitoring it on input nodes its not used anywhere else
  neuron.outputNumber = 0 --this is just for monitoring it on output nodes its not used anywhere else
	neuron.incoming = {} --data from previous thingis
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
    if  randN > geneActivationChance then
      tempConnectionGene.status = true
      print("One or more of connection genes are enabled")
    else
      tempConnectionGene.status = false
      print("One or more of connection genes are disabled")
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
    print("old neuron value".. genome.network[i].value)
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
         print("Associated gene in"..genome.genes[k].input.value)
         print("Associated WEIGHT "..genome.genes[k].weight)
         print("Associated gene out"..genome.genes[k].out.value)
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
     -- print("sum: "..sum)
      end
      activation = sigmoid(sum)
      --print("activated")
    --  print("activation value: "..activation)
      if #tempWout~=0 then
        genome.network[i].value = activation
      --replace all values of #tempwout.out with new activation value
      for n = 1, #tempWout do
        tempWout[n].out.value = activation
      end
        end

      --replace in actual gene as well (both ins and outs)
      --REPLACE NODE
    print("new neuron value".. genome.network[i].value)
    --print("new genome i/o stats".. genome.network[i].inStatus)
  end
  obtainOutputs(genome)
end

function createStartingPopulation(number)
  genomesCreated = {}
  for i = 1, number do
    testPop = createNewGenome()
    BuildNetwork(testPop)
    if nodeMutationChance > math.random() then
      mutateNodeGene(testPop)
    end
    if geneMutationChance > math.random() then
      mutateConnectionGene(testPop)
    end
    --evaluateNetwork(testPop)
    table.insert(genomesCreated,testPop)
    print("population no: "..i)
  end
  return genomesCreated
  end

--local initialPop = createStartingPopulation(initialPopulationSize)
--print("population size"..#initialPop)
--the onl neurons it should be accepting are output neurons
testPop = createNewGenome()
BuildNetwork(testPop)
evaluateGenome(testPop)
print("inStatus "..testPop.network[9].inStatus)


c = 0
print("c stats "..#testPop.network)
for i = 1, #testPop.network do
    if testPop.network[i].inStatus == 1 then
      c = c + 1
      end
end
print("c stats "..c)
function runGeneration(NumberOfGenerations)
  for i = 1, NumberOfGenerations do
    for j = 1, #initialPop do
      evaluateNetwork(initialPop[j])

    end
  end
end

mutateConnectionGene(testPop)
mutateNodeGene(testPop)
evaluateGenome(testPop)
mutateNodeGene(testPop)
mutateNodeGene(testPop)
mutateNodeGene(testPop)
mutateNodeGene(testPop)
mutateNodeGene(testPop)
evaluateGenome(testPop)
mutateNodeGene(testPop)
mutateNodeGene(testPop)
mutateNodeGene(testPop)
mutateConnectionGene(testPop)
mutateConnectionGene(testPop)
mutateConnectionGene(testPop)
evaluateGenome(testPop)
mutateNodeGene(testPop)
mutateNodeGene(testPop)
mutateNodeGene(testPop)
mutateConnectionGene(testPop)
mutateConnectionGene(testPop)
mutateConnectionGene(testPop)
evaluateGenome(testPop)
mutateNodeGene(testPop)
mutateNodeGene(testPop)
mutateNodeGene(testPop)
mutateConnectionGene(testPop)
mutateConnectionGene(testPop)
mutateConnectionGene(testPop)
evaluateGenome(testPop)
mutateNodeGene(testPop)
mutateNodeGene(testPop)
mutateNodeGene(testPop)
mutateConnectionGene(testPop)
mutateConnectionGene(testPop)
mutateConnectionGene(testPop)
evaluateGenome(testPop)
--print(sigmoid(0))
