 local TestInputs = {}
TestInputs[1] = 12
TestInputs[2] = 39
TestInputs[3] = 34
TestInputs[3] = 629
TestInputs[3] = 629

local TestOutputs = {}
TestOutputs[1] = 0
TestOutputs[2] = 01
TestOutputs[3] = 03
TestOutputs[4] = 04

nodeMutationChance = 0.5  --in use
geneMutationChance = 0.2  --in use
geneActivationChance = 0.5
selectionChance = 0.3 --in use
initialPopulationSize = 1000 --in use
speciesDistance = 0.23 --in use
crossOverChance = 0.75


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
print("max val in connection gene"..max)
node1 = genome.genes[math.random(1,#genome.genes)]
node2 = genome.genes[math.random(1,#genome.genes)]
state1 = 0
state2 = 0
if node1.input ~= node2.out and node1.input.inStatus ~= 0 then
state1 = 1
end
if node1.out~=node2.input and node1.out.inStatus ~= 1 then
state2 = 1
end
if state1 == 1 and state2 == 1 then
connectionGeneMutate1 = connectionGene()
connectionGeneMutate1.input = node1.input
connectionGeneMutate1.innovation = max + 1
table.insert(connectionGeneMutate1.input.weightIndex,connectionGeneMutate1.innovation)
connectionGeneMutate1.out = node2.input
table.insert(connectionGeneMutate1.out.weightIndex,connectionGeneMutate1.innovation)
print("a gene has been mutated by connection")
table.insert(genome,connectionGeneMutate1)
end
return genome
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
connectionGeneMutate1 = connectionGene()
connectionGeneMutate1.input = connectionGeneTemp.input --attatch neuron
connectionGeneMutate1.innovation = maxInnovation + 1
table.insert(connectionGeneMutate1.input.weightIndex,connectionGeneMutate1.innovation)
tempNeuron = newNeuron()
tempNeuron.value = math.random()
table.insert(tempNeuron.weightIndex,connectionGeneMutate1.innovation)
connectionGeneMutate1.out = tempNeuron
table.insert(genome.network,tempNeuron)
table.insert(genome.genes,connectionGeneMutate1)
connectionGeneMutate2 = connectionGene()
connectionGeneMutate2.input = tempNeuron
connectionGeneMutate2.innovation = connectionGeneMutate1.innovation + 1
table.insert(tempNeuron.weightIndex,connectionGeneMutate2.innovation)
connectionGeneMutate2.out = connectionGeneTemp.out
table.insert(connectionGeneMutate2.out.weightIndex,connectionGeneMutate2.innovation)
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
	for i=1,#TestInputs do
    tempN = newNeuron()
    tempN.value = TestInputs[i]
    tempN.inStatus = 1 --signifies input neuron
    tempN.inputNumber = i
    --create gene and link innovation number
    tempConnectionGene = connectionGene()
    tempConnectionGene.input = tempN
    tempConnectionGene.innovation = innovationNumber
    if math.random() < geneActivationChance then
      tempConnectionGene.status = true
      print("One or more of connection genes are enabled")
    else
      tempConnectionGene.status = false
      print("One or more of connection genes are disabled")
      end
    --link neuron
    table.insert(tempN.weightIndex,innovationNumber)
    table.insert(genome.genes,tempConnectionGene)
    --add to genomee network
    table.insert(genome.network,tempN)
    innovationNumber = innovationNumber + 1
	end

	--create outputs
	for o=1,#TestOutputs do
    tempO = newNeuron()
    tempO.value =  TestOutputs[o]
    tempO.inStatus = 0 --outputneuron
    tempN.outputNumber = o
    --loop through genes and link all "outs" to this output
    for i = 1, #genome.genes do
      --link neuron to genes through innovation number in genes
      table.insert(tempO.weightIndex,genome.genes[i].innovation)
      genome.genes[i].out = tempO
      end
    --add neuron to network
    table.insert(genome.network,tempO)
	end
  --return genome
end

function updateInputs(genome)
  inputCount = 1
  for i = 1, #genome.network do
    if genome.network[i].inStatus==1 and inputCount <= #TestInputs then
      --update input var
      genome.network[i].value = TestInputs[inputCount]
      inputCount = inputCount + 1
    end
  end
end

--for testing purposes
  function obtainOutputs(genome)
    OutCount = 1
  for i = 1, #genome.network do
    if genome.network[i].inStatus==0 and OutCount <= #TestOutputs then
      --update input var
      print("New Output"..genome.network[i].value)
      OutCount = OutCount + 1
      end
    end
    end

function evaluateNetwork(genome)
  --obtain/update input neurons (input neurons have a state of 1)
  updateInputs(genome)
  --obtain all connections to a node and shit out output
  for i = 1, #genome.network do
    tempW = {} --all associated (according to innovation numbers)
    tempWout = {} --just the outs
    oldValue = genome.network[i] --old neuron data
    --update node AND gene
    --obtain node value
    --TempNode = genome.network[i]
    print("genome network neuron old::  " ..genome.network[i].value)
    print("weight index values " ..#genome.network[i].weightIndex) --how many innovation numbers are in here
    --obtain weights through index
    for j = 1, #genome.network[i].weightIndex do
      --grab innovation number and loop through network weights to check if  there is a match
      tempWeightIndx = genome.network[i].weightIndex[j]
      print("innovation number to look for"..tempWeightIndx)
      for k = 1, #genome.genes do
        --if innovation number of neuron is equal to innovation number of connection gene
        if tempWeightIndx == genome.genes[k].innovation then
          --clear a table: for k,v in pairs(tab) do tab[k]=nil end
          --add to tempW
          print("Associated gene in"..genome.genes[k].input.value)
          print("Associated gene in"..genome.genes[k].out.value)
          --add to tempW(all ins and outs connection genes)
          table.insert(tempW,genome.genes[k])--this will now hold all associated connection genes
          end
        end
      end
      --after finding all associated genes, filter ins and outs(if gene out == genome.network) and store temporarily
      print("number of all found matching genes"..#tempW)
      for b = 1, #tempW do
        print("gene status: "..tostring(tempW[b].status))
        print("gene out: "..tempW[b].out.value)
        print("network out: "..genome.network[i].value)
        if tempW[b].out== genome.network[i] and tempW[b].status == true then
          print("matched actual neurons")
          --now grab gene and store it
          table.insert(tempWout, tempW[b])
        end
        if tempW[b].out== genome.network[i] and tempW[b].status == false then
          print("right output but set to inactive so wont be processed")
        end
    end
    print("found refined connection genes : "..#tempWout)
    --loop through all cases of tempWout, obtain ins, multiply by weights get new out value and REPLACE the gene with the new tempWouts
    sum = 0
    activation = 0
    for m = 1, #tempWout do
      print("input value from filtered "..tempWout[m].input.value)
      sum = sum + (tempWout[m].input.value * tempWout[m].weight)
    end
    if sum~=0 then
      print("activated")
      activation = sigmoid(sum)
      print("activation: "..activation)
      genome.network[i].value = activation
      --replace all values of #tempwout.out with new activation value
      for n = 1, #tempWout do
        tempWout[n].out.value = activation
      end
      --replace in actual gene as well (both ins and outs)
      --REPLACE NODE
    end
    print("new genome value".. genome.network[i].value)
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

initialPop = createStartingPopulation(initialPopulationSize)
print("population size"..#initialPop)

