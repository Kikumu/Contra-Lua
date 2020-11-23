 playerlocationYaxis = memory.readbyte(0x031A)
  playerlocationXaxis = memory.readbyte(0x0334)
  playerHealth = memory.readbyte(0x0032)
  enemylocationXaxis1 = memory.readbyte(0x0343)
  enemylocationXaxis2 = memory.readbyte(0x0344)
  screenMov1 = memory.readbyte(0x0065)
  screenMov2 = memory.readbyte(0x00FD)



local TestInputs = {}
TestInputs[1] = playerlocationYaxis
TestInputs[2] = playerlocationXaxis
TestInputs[3] = playerHealth
TestInputs[4] = enemylocationXaxis1
TestInputs[5] = enemylocationXaxis2
TestInputs[6] = screenMov1
TestInputs[7] = screenMov2

local TestOutputs = {}
TestOutputs[1] = 0
TestOutputs[2] = 0
TestOutputs[3] = 0
TestOutputs[4] = 0
TestOutputs[5] = 0
TestOutputs[6] = 0
TestOutputs[7] = 0
TestOutputs[8] = 0

ButtonNames = {
  "A",
  "B",
  "X",
  "Y",
  "Up",
  "Down",
  "Left",
  "Right",
}


function clearJoypad()
	controller = {}
	for b = 1,#ButtonNames do
		controller["P1 " .. ButtonNames[b]] = false
	end
	joypad.set(controller)
end


--loop through output
--check/look for out put in network
--activate/deactivate output value
function moveAnimation(genome)
  for i = 1, #TestOutputs do
    for k = 1, #genome.network do
      if genome.network[k] > 0 and genome.network[k].inStatus == 0 and genome.network[k].outputNumber == i then
        highRoller = {}
        highRoller["P1 "..ButtonNames[i]] = true
        joypad.set(highRoller)
      end
      if genome.network[k] < 0 and genome.network[k].inStatus == 0 and genome.network[k].outputNumber == i then
        highRoller = {}
        highRoller["P1 "..ButtonNames[i]] = false
        joypad.set(highRoller)
      end
    end
    end
end


nodeMutationChance = 0.5  --in use
geneMutationChance = 0.2  --in use
geneActivationChance = 0.5
selectionChance = 0.3 --in use
initialPopulationSize = 1000
speciesDistance = 0.23
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
--[[there is one thing that bothers me...i have physical locations of every connection gene but no
actual place for the "data points" will correct this but first..some cod!!! DONE
]]
--[[

TO DO:
1) Take it out for a test spin to test out functions created --DONE
2) Test propagation --Done
3) Neurons holder for each network created --DONE
4) Find a way to keep track of the fitness of each genome --done
5) Create function to organise a group of genes into different speciesPool
6) Rework how the network breeds --done
7) Fix how genes mutate --done


]]

local neurons = {}
neurons[1] = 20

--activation function
function sigmoid(x)
    return 1.0 / (1.0 + math.exp(-x))
end



--how many inputs and outputs? 3 inputs and one output to start with
--all it does for now is just create connections to the main input node
function connectionGene()
--this is like the weights of the network instructing on what connects with what
local gene = {}

gene.input = 0
gene.out = 0
gene.weight = math.random()
gene.status = true
gene.innovation = 0 --ancestry monitor

return gene

end

--a((i*w) + (i*w) + (i*w))


--tested(needs reviewing) --TESTED AND REWORKED
function mutateConnectionGene(genome)
--take 2 random nodes and adds a connection between them if none are available
max = retunMaxInnovation(genome)
print("max val in connection gene"..max)
node1 = genome.genes[math.random(1,#genome.genes)]
node2 = genome.genes[math.random(1,#genome.genes)]
--print ('Node 1 '..node1.weight)
--print ('Node 2 '..node2.weight)
state1 = 0
state2 = 0

--check if there is a connection.(in in or out) if no connection, add one and make 2 new connection genes
--do not initiate a connection if the input node is an output
if node1.input ~= node2.out and node1.input.inStatus ~= 0 then
state1 = 1
end
--do not initiate a connection if the output node is an input
if node1.out~=node2.input and node1.out.inStatus ~= 1 then
state2 = 1
end

--
if state1 == 1 and state2 == 1 then
connectionGeneMutate1 = connectionGene()
connectionGeneMutate1.input = node1.input
--add connection to neuron
connectionGeneMutate1.innovation = max + 1
table.insert(connectionGeneMutate1.input.weightIndex,connectionGeneMutate1.innovation)
connectionGeneMutate1.out = node2.input
--add connection to neuron at gene out
table.insert(connectionGeneMutate1.out.weightIndex,connectionGeneMutate1.innovation)
print("a gene has been mutated by connection")
table.insert(genome,connectionGeneMutate1)
return genome
end

--if state2 == 1 then

--end

print('state1: '..state1)
print('state2: '..state2)
return gene
end

--TESTED AND REWORED
function retunMaxInnovation(gene)
maxInnovation = gene.genes[1].innovation
print("Max Innovation: "..maxInnovation)
for i = 1, #gene.genes do
  if gene.genes[i].innovation > maxInnovation then
  print("looping innovation in genes.."..gene.genes[i].innovation)
  maxInnovation = gene.genes[i].innovation
  end
end

return maxInnovation
end

--FULLY REWORKED
function mutateNodeGene(genome)
--pick a random node/neuron/connection gene
--max innovation
tstval = #genome.genes
maxInnovation = retunMaxInnovation(genome)
print("max innovation"..maxInnovation)
print("Initial length of gene: "..#genome.genes)
genepos = math.random(1,#genome.genes)

--i dont want a connection gene going out of an output

print("selected random pos for mutation point"..genepos)
connectionGeneTemp = genome.genes[genepos]
print("gene connection gene status: "..tostring(genome.genes[genepos].status))
print(" neuron point value to be mutated input".. connectionGeneTemp.input.value)
print(" neuron point value to be mutated output".. connectionGeneTemp.out.value)
--disable gene
connectionGeneTemp.status = false
print("Temp connection gene status: "..tostring(connectionGeneTemp.status))
print("gene connection gene status: "..tostring(genome.genes[genepos].status))
--put it back
--genome.genes[genepos] = connectionGeneTemp (no need i guess)...
--print("Length check after putting it back: "..#genome.genes)
--print("Original innovation"..node1.innovation)
--node one in to new connection to new neuron

--crate new connection
connectionGeneMutate1 = connectionGene()
connectionGeneMutate1.input = connectionGeneTemp.input --attatch neuron
print("new connection gene input: "..connectionGeneMutate1.input.value)
connectionGeneMutate1.innovation = maxInnovation + 1
--add gene innovation number to neuron of input
table.insert(connectionGeneMutate1.input.weightIndex,connectionGeneMutate1.innovation)
print("number of values in this neurons index"..#connectionGeneMutate1.input.weightIndex)

print("new connection gene innovation: "..connectionGeneMutate1.innovation)
print("number of neurons in my gene"..#genome.network)
--create a neuron and put into network
tempNeuron = newNeuron()
tempNeuron.value = math.random()
--add weight index
table.insert(tempNeuron.weightIndex,connectionGeneMutate1.innovation)
--tempNeuron.weightIndex = connectionGeneMutate1.innovation
connectionGeneMutate1.out = tempNeuron
--add neuron to network
table.insert(genome.network,tempNeuron)
print("new connection gene output"..connectionGeneMutate1.out.value)
table.insert(genome.genes,connectionGeneMutate1)
print("connection gene in: "..genome.genes[#genome.genes].input.value)
print("connection gene out: "..genome.genes[#genome.genes].out.value)
print("Gene length after adding a new connection and gene: "..#genome.genes)

connectionGeneMutate2 = connectionGene()
connectionGeneMutate2.input = tempNeuron
print("new 2nd connection gene input: "..connectionGeneMutate2.input.value)
connectionGeneMutate2.innovation = connectionGeneMutate1.innovation + 1
--link neuron 2 weight by adding innovation number to neuron
table.insert(tempNeuron.weightIndex,connectionGeneMutate2.innovation)

print("new 2nd connection gene innovation: "..connectionGeneMutate2.innovation)
connectionGeneMutate2.out = connectionGeneTemp.out
--add connection gene innovation no. to out neuron of connection gene
table.insert(connectionGeneMutate2.out.weightIndex,connectionGeneMutate2.innovation)
print("new 2nd connection gene out: "..connectionGeneMutate2.out.value)
--add gene to network
table.insert(genome.genes,connectionGeneMutate2)

print("Gene length after adding a final connection and gene: "..#genome.genes)
--return gene
end



--for crossover purposes(Goal is to add to g1)tested--TESTED AND REWORKED
function matchingGenes(genome1,genome2)
matchedGenes ={}
neuronsExtracted = {}
g1Length = #genome1.genes
g2Length = #genome2.genes

--check innovations. if innovations match, add to bucket
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
--return matched genes containing the same innovation number(some may be active or inactive)
return matchedGenes, neuronsExtracted
end


--(takes disjointed genes from g2), goal is to add to G1
--the genes in the middle
function disjointGenes(genome1,genome2)

--pick max innovation number from g1 asnd g2
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


--confident(tested)
function excessGenes(genome1,genome2)
ExcessGenesTable = {}
neuronsExtracted = {}
--find excess genes in g2(which are located in g1)
maxInnovationg1 =retunMaxInnovation(genome1)
print("max innovation in g1 in excess genes monitor"..maxInnovationg1)
maxInnovationg2 = retunMaxInnovation(genome2)
print("max innovation in g2 in excess genes monitor"..maxInnovationg2)
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
print("in selection sort,d to sortNumber "..#DisjointedGenesArr)
print("in selection sort,e to sortNumber "..#ExcessGenesArr)
print("in selection sort,m to sortNumber "..#MatchingGenesArr)

--table.insert(genesToSort,DisjointedGenesArr)
--table.insert(genesToSort,ExcessGenesArr)
--table.insert(genesToSort,MatchingGenesArr)
accumilateGenesForSorting(DisjointedGenesArr,genesToSort)
accumilateGenesForSorting(ExcessGenesArr,genesToSort)
accumilateGenesForSorting(MatchingGenesArr,genesToSort)


print("in selection sort,genes to sortNumber "..#genesToSort)
for i = 1, #genesToSort do
minInnovationIndex = i
for j = i + 1, #genesToSort do
print("selection no: "..j)
print("innovation in sort..........."..genesToSort[j].innovation)
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
  --total = #n1+#n2+#n3
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
--childGenome= {}
childlength = 0
--gene1 = genome.genes[math.random(1,#genome)]
--gene2 = genome.genes[math.random(1,#genome)]
m,n1 = matchingGenes(genome1,genome2)
d,n2 = disjointGenes(genome1,genome2)
e,n3 = excessGenes(genome1,genome2)

nT = combineNeurons(n1,n2,n3)
nT1 = removeDuplicateNeurons(nT)
print("Total number of neurons after crossover"..#nT)
print("Total number of clean neurons after crossover"..#nT1)
genome = selectionSort(d, e, m)
genome.network = nT1
return genome

end



function species(genes)

speciesPool = {}

return speciesPool
end

--neurons should match with innovation number
--this is my new 'newGenome' function
function createNewGenome()
  genomeCluster = {} --holds the gene information network (general info)
  genomeCluster.genes = {} --weight info (connection genes info)
  genomeCluster.fitness = 0
  genomeCluster.network = {} --holds neurons
  --genomeCluster.weightIndex = {} --holds index of a node that its connected with per neuron array

  return genomeCluster
  end
--testPropagate = newGenome()
--try and forward propagate testPropagate
--find all connected lines and

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

--you only need to do this ONCE
 function BuildNetwork(genome)
	print("network size initial value"..#genome.network)
  --5 inputs and 4 outputs
  innovationNumber = 1
  --add neurons(create or direct inputs)
	for i=1,#TestInputs do
    tempN = newNeuron()
    tempN.value = TestInputs[i]
    tempN.inStatus = 1 --signifies input neuron
    tempN.inputNumber = i
    --create gene and link innovation number
    tempConnectionGene = connectionGene()
    tempConnectionGene.input = tempN
    tempConnectionGene.innovation = innovationNumber
    --link neuron
    table.insert(tempN.weightIndex,innovationNumber)
    --remember--output of this gene is nil...for now..
    --insert connection gene back 2 de genome
    table.insert(genome.genes,tempConnectionGene)
    --table.insert(network.neurons,tempN)
    --add to genomee network
    table.insert(genome.network,tempN)
    innovationNumber = innovationNumber + 1
   -- print("neuron input value: "..network.neurons[i].value)
	end

  --if #genome.network > #inputs + #outputs then
  --remember weight index in neuron links to innovation number of gene
   -- end


	--create outputs
	for o=1,#TestOutputs do
    tempO = newNeuron()
    tempO.value =  TestOutputs[o]
    tempO.inStatus = 0 --outputneuron
    tempO.outputNumber = o
    --loop through genes and link all "outs" to this output
    for i = 1, #genome.genes do

      --link neuron to genes through innovation number in genes
      table.insert(tempO.weightIndex,genome.genes[i].innovation)
      genome.genes[i].out = tempO
      --print("gene input xd: "..genome.genes[i].input.value)
      --print("gene output xd: "..genome.genes[i].out.value)
      end
    --add neuron to network
    table.insert(genome.network,tempO)
    -- print("neuron output value: "..network.neurons[o].value)
	end


  print("network size after adding inputs and outputs value"..#genome.network)
  return genome
end

function updateInputs(genome)
  --loop through input
  --loop through network and obtain matching input
  --update accordingly
  for i = 1, #TestInputs do
    for j = 1, #genome.network do
      if genome.network[j].inputNumber == i then
        genome.network[j].input.value = TestInputs[i]
      end
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
    testPop = BuildNetwork(testPop)
    if nodeMutationChance > math.random() then
      mutateNodeGene(testPop)
    end
    if geneMutationChance > math.random() then
      mutateConnectionGene(testPop)
    end
    evaluateNetwork(testPop)
    table.insert(genomesCreated,testPop)
  end
  return genomesCreated
  end

initialPop = createStartingPopulation(initialPopulationSize)
print("population size"..#initialPop)
--for i = 1, #initialPop do
  --pick random gene as a mascot
--  mascot = initialPop[math.random(1,#initialPop)]
--  speciation(initialPop[i],mascot)
--  end
