mutationChance = 0.5
geneActivationChance = 0.5
selectionChance = 0.3 --this value is used to choose between 2 connection genes of the same innovation
--species classification variables

disjointment = 0.2 --c1
excessGenes = 0.3 --c2
weightImportance = 0.1 --c3
MaxNodes = 11

--E is number of disjointed connection genes
--D is number of excess connection genes
--W is weight value
--N number of connection genes

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
6) Rework how the network breeds
7) Some polish I guess


]]

local neurons = {}
neurons[1] = 20


local inputs = {}
inputs[1] = 12
inputs[2] = 39
inputs[3] = 69

local outputs = {}
outputs[1] = 0

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


--tested(needs reviewing)
function mutateConnectionGene(gene)
--take 2 random nodes and adds a connection between them if none are available
max = retunMaxInnovation(gene)
node1 = gene[math.random(1,#gene)]
node2 = gene[math.random(1,#gene)]
--print ('Node 1 '..node1.weight)
--print ('Node 2 '..node2.weight)
state1 = 0
state2 = 0

--check if there is a connection.(in in or out) if no connection, add one and make 2 new connection genes
if node1.input ~= node2.out then
state1 = 1
end

if node1.out~=node2.input then
state2 = 1
end

--
if state1 == 1 and state2 == 1 then
connectionGeneMutate1 = connectionGene()
connectionGeneMutate1.input = node1.input
connectionGeneMutate1.innovation = max + 1
connectionGeneMutate1.out = node2.input
--connectionGeneMutate1.output = a(i*w)?
print("a gene has been mutated by connection")
connectionGeneMutate1.out= node2.input
table.insert(gene,connectionGeneMutate1)

return gene
end

--if state2 == 1 then

--end

print('state1: '..state1)
print('state2: '..state2)
return gene
end


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


function mutateNodeGene(genome)
--pick a random node/neuron/connection gene
--max innovation
maxInnovation = retunMaxInnovation(genome)
print("max innovation"..maxInnovation)
print("Initial length of gene: "..#genome.genes)
genepos = math.random(1,#gene.genes)
print("selected random pos for mutation point"..genepos)
node1 = genome.genes[genepos]
print(" point to be mutated input".. node1.input)
print(" point to be mutated output".. node1.out.value)
--disable gene
node1.status = false
--put it back

genome.genes[genepos] = node1
print("Length check after putting it back: "..#genome.genes)
--print("Original innovation"..node1.innovation)
--node one in to new connection to new neuron

connectionGeneMutate1 = connectionGene()

connectionGeneMutate1.input = node1.input
print("new connection gene input: "..connectionGeneMutate1.input)
connectionGeneMutate1.innovation = maxInnovation + 1
print("new connection gene innovation: "..connectionGeneMutate1.innovation)
print("number of neurons in my gene"..#genome.genes.network)
--create a neuron and put into network
tempNeuron = newNeuron()
tempNeuron.value = math.random(60,200)
tempNeuron.weightIndex = connectionGeneMutate1.innovation
table.insert(genome.genes.network,tempNeuron)
connectionGeneMutate1.out = tempNeuron
print("new connection gene output"..connectionGeneMutate1.out.value)
table.insert(genome.genes,connectionGeneMutate1)

print("Gene length after adding a new connection and gene: "..#genome.genes)

connectionGeneMutate2 = connectionGene()
connectionGeneMutate2.input = genome.genes.network[#genome.genes.network]
print("new 2nd connection gene input: "..connectionGeneMutate2.input.value)
connectionGeneMutate2.innovation = connectionGeneMutate1.innovation + 1
print("new 2nd connection gene innovation: "..connectionGeneMutate2.innovation)
connectionGeneMutate2.out = node1.out
print("new 2nd connection gene out: "..connectionGeneMutate2.out)
table.insert(genome.genes,connectionGeneMutate2)

print("Gene length after adding a final connection and gene: "..#gene)

return gene
end



--for crossover purposes(Goal is to add to g1)tested
function matchingGenes(g1,g2)
matchedGenes ={}
g1Length = #g1
g2Length = #g2

for i = 1, g1Length do
for j = 1, g2Length do
if g1[i].innovation == g2[j].innovation then
if selectionChance > math.random() then
table.insert(matchedGenes,g1[i])
else
table.insert(matchedGenes,g2[j])
end
end
end
end

return matchedGenes
end


--(takes disjointed genes from g2), goal is to add to G1
function disjointGenes(g1,g2)

--pick max innovation number from g1 asnd g2
maxInnovation = 0
disjointedGenes = {}
g2MaxInnovation = retunMaxInnovation(g2)
g1MaxInnovation = retunMaxInnovation(g1)
g2Length = #g2
g1Length = #g1

found = 0
for i = 1, #g2 do
  for j = 1, #g1 do
    if g1[j].innovation < g2MaxInnovation then
    if g2[i].innovation == g1[j].innovation then
      found = 1
    end
  end
if found == 0 then
table.insert(disjointedGenes,g2[i])
end
    end

end
return disjointedGenes

end


--confident(tested)
function excessGenes(g1,g2)
ExcessGenesTable = {}
--find excess genes in g2(which are located in g1)
maxInnovationg1 = g1[#g1].innovation
maxInnovationg2 = g2[#g2].innovation
g1Length = #g1
g2Length = #g2

if maxInnovationg1 > maxInnovationg2 then
for i = 1, g1Length do
if g1[i].innovation > maxInnovationg2 then
table.insert(ExcessGenesTable, g1[i])
end
end
return ExcessGenesTable
end

if maxInnovationg2 > maxInnovationg1 then
for i = 1, g2Length do
if g2[i].innovation > maxInnovationg1 then
table.insert(ExcessGenesTable, g2[i])
end
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

return genesToSort
end


--breed
function crossover(genome) --genome in the sense that you are passing to this function a set of genomes
--childGenome= {}
childlength = 0
gene1 = genome[math.random(1,#genome)]
gene2 = genome[math.random(1,#genome)]
m = matchingGenes(gene1,gene2)
d = disjointGenes(gene1,gene2)
e = excessGenes(gene1,gene2)
res = selectionSort(d, e, m)
return res

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
  neuron.inStatus = {} --0 if an input 1 if an output (matches with node array index
	neuron.incoming = {} --data from previous thingis
	neuron.value = 0.0 -- current neuron value
	return neuron
end

--gene.input = 0
--gene.out = 0
--gene.weight = math.random()
--gene.status = true
--gene.innovation = 0 --ancestry monitor

 function BuildNetwork(genome)
	print("network size initial value"..#genome.network)

  innovationNumber = 1
  --add neurons(create or direct inputs)
	for i=1,#inputs do
    tempN = newNeuron()
    tempN.value = math.random(200,300)

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
	for o=1,#outputs do
    tempO = newNeuron()
    tempO.value =  math.random(500,800)
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


function evaluateNetwork(genome)
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
      for p = 1, #tempWout do
        for q = 1, #genome.genes do
          --if same innovation number and same output
          if genome.genes[q].innovation == tempWout[p].innovation and genome.genes[q].out == oldValue then
            print("previous "..genome.genes[q].out.value)
            genome.genes[q]= tempWout[p]
            print("current "..genome.genes[q].out.value)
          end
          --if same innovation number and same input
          if genome.genes[q].innovation == tempWout[p].innovation and genome.genes[q].input == oldValue then
            print("previous "..genome.genes[q].input.value)
            genome.genes[q]= tempWout[p]
            print("current "..genome.genes[q].input.value)
          end

          end
        end
      --REPLACE NODE
    end
    print("new genome value".. genome.network[i].value)
  end
  end
 testPropagate = createNewGenome()

 testPropagate = BuildNetwork(testPropagate)

 --testPropagate = mutateNodeGene(testPropagate)

 f_to_pay_respects = evaluateNetwork(testPropagate)
 --print("gene value"..testPropagate.network[2].value)
 --print("gene out value"..testPropagate.genes[2].out.value)
 --I want to find out how many weight values are attatched to me



 --WEIGHT INDEXES AREE FOUND IN NEURONS U DUMBASS

--  function newNeuron()
--	local neuron = {}
--  neuron.weightIndex = {} --keeps track of attatched weights (matches with node array index
--  neuron.inStatus = {} --0 if an input 1 if an output (matches with node array index
--	neuron.incoming = {} --data from previous thingis
--	neuron.value = 0.0 -- current neuron value
--	return neuron
--end
