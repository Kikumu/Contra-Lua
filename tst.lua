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



--FIND A WAY TO RESTRUCTURE THE WAY OU START CREATING YOUR GENES. FROM CURRENT SPECULATION, THE POPULATION CREATION WORKS ALONGSIDE ALL THESE FACTORS TO CREATE DIFFERENT SPECIES (HENCE SPECIES POOL WITHIN THE POPULATION)
function newGenome() --just creates starting pop, a single genome

innovationGene = 0
newGeneOut = {}

for i = 1,#inputs do
parentgene = connectionGene()
parentgene.input = inputs[i]
parentgene.innovation = innovationGene
parentgene.out = outputs[1]
table.insert(newGeneOut,i,parentgene)
innovationGene = innovationGene + 1
end

return newGeneOut

end


--tested
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
maxInnovation = gene[1].innovation
for i = 1, #gene do
  if gene[i].innovation > maxInnovation then
  maxInnovation = gene[i].innovation
  end
end

return maxInnovation
end


function mutateNodeGene(gene)
--pick a random node/neuron/connection gene
--max innovation
maxInnovation = retunMaxInnovation(gene)
print("max innovation"..maxInnovation)
print("Initial length of gene: "..#gene)
genepos = math.random(1,#gene)
print("selected random pos for mutation point"..genepos)
node1 = gene[genepos]
print(" point to be mutated input".. node1.input)
print(" point to be mutated output".. node1.out)
--disable gene
node1.status = false

--put it back

gene[genepos] = node1
print("Length check after putting it back: "..#gene)
--print("Original innovation"..node1.innovation)
--node one in to new connection to new neuron

connectionGeneMutate1 = connectionGene()

connectionGeneMutate1.input = node1.input
print("new connection gene input: "..connectionGeneMutate1.input)
connectionGeneMutate1.innovation = maxInnovation + 1
print("new connection gene innovation: "..connectionGeneMutate1.innovation)
table.insert(neurons,math.random(60,200))
connectionGeneMutate1.out = neurons[#neurons]
print("new connection gene output"..connectionGeneMutate1.out)
table.insert(gene,connectionGeneMutate1)

print("Gene length after adding a new connection and gene: "..#gene)

connectionGeneMutate2 = connectionGene()
connectionGeneMutate2.input = neurons[#neurons]
print("new 2nd connection gene input: "..connectionGeneMutate2.input)
connectionGeneMutate2.innovation = connectionGeneMutate1.innovation + 1
print("new 2nd connection gene innovation: "..connectionGeneMutate2.innovation)
connectionGeneMutate2.out = node1.out
print("new 2nd connection gene out: "..connectionGeneMutate2.out)
table.insert(gene,connectionGeneMutate2)

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


--2 mutations therefore the network grew by a size of 2 (total 7)


xArr = {}
yArr = {}

dpoints = #inputs + #outputs + #neurons
--dInputs = #inputs
--dOutput = #outputs
--dNeurons= #neurons



math.randomseed(os.time())
for i = 1, dpoints do
  table.insert(xArr,math.random(1,50)* 10)
  --print("Value: "..xArr[i])
end

for i = 1, dpoints do
  table.insert(yArr,math.random(1,100) * 3)
  --print("Value: "..yArr[i])
end
--xpos = math.random(20,100)
--draw 5 neurons in 5 random places
function drawLine(x1,y1,x2,y2)

  end

--love.graphics.line( x1, y1, x2, y2,)
xInputArr = {}
yInputArr = {}
--for i = 1, dInputs do
 -- table.insert(xInputArr,math.random(20,100)* 5 )
--  table.insert(yInputArr,math.random(20,100) * 3)
--  end

inputsX = {}
inputsY = {}
inputsObtained = {}
inputMonitor = 1


--search for neurons and assign a "draw"

nv = #inputs + #neurons

countNodes = 1
for i = 1,#inputs do
  table.insert(inputsObtained,inputs[i])
  table.insert(inputsX,xArr[countNodes])
  table.insert(inputsY,xArr[countNodes])
  inputsX[countNodes] = xArr[countNodes]
  inputsY[countNodes] = yArr[countNodes]
  countNodes = countNodes + 1
end
for i = 1, #neurons do
  table.insert(inputsObtained,neurons[i])
  table.insert(inputsX,xArr[countNodes])
  table.insert(inputsY,xArr[countNodes])
  inputsX[countNodes] = xArr[countNodes]
  inputsY[countNodes] = yArr[countNodes]
  countNodes = countNodes + 1
end
for i = 1, #outputs do
  table.insert(inputsObtained,outputs[i])
  table.insert(inputsX,xArr[countNodes])
  table.insert(inputsY,xArr[countNodes])
  inputsX[countNodes] = xArr[countNodes]
  inputsY[countNodes] = yArr[countNodes]
  countNodes = countNodes + 1
end


--love.graphics.line( x1, y1, x2, y2)

print("number of values " ..countNodes)
function love.draw()
  for i = 1, nv do
    love.graphics.setColor(1, 0.5, 1)
    love.graphics.circle("fill", xArr[i], yArr[i], 10, 100) -- Draw white circle with 100 segments.
    love.graphics.setColor(0, 1, 0.3, 1)
    love.graphics.print(""..inputsObtained[i], inputsX[i], inputsY[i])
    --love.graphics.line( x1, y1, x2, y2)
  end
end




--print(""..neurons[2])
--print("number of yval " ..y[5].input)
geneCluster = {}
x = newGenome()
mutateNodeGene(x)
mutateNodeGene(x)

y = newGenome()
mutateNodeGene(y)
mutateConnectionGene(y)

--
table.insert(geneCluster,x)
table.insert(geneCluster,y)

print("gene cluster "..#geneCluster)
z = crossover(geneCluster)

for i=1,#x do
print("in"..x[i].input)
print("out"..x[i].out)
end
print("gene 1 "..#x)
--print("gene 2 "..#y)
for i=1,#y do
print("in"..y[i].input)
print("out"..y[i].out)
end
print("gene 2 "..#y)
print("gene 3 "..#z)
--a = matchingGenes(x,y)
--a = excessGenes(x,y)
--a = disjointGenes(x,y)
print("number of genes in child"..#z)
for i=1,#z do
print("in"..z[i].input)
print("out"..z[i].out)
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
  neuron.weightIndex = {} --keeps track of attatched weights (matches with node array index
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

 function createNewNetwork(genome)
  local network = {}
	network.neurons = {} --previous inputs (and weights mabe)?
  network.weights = {}
  network.weightIndex = {}
	print("network size initial value"..#network.neurons)

  weightTracker = 1
  --add neurons(create or direct inputs)
	for i=1,#inputs do
    tempN = newNeuron()
    tempN.value = math.random(200,300)
    table.insert(tempN.weightIndex,weightTracker) --stores tracker to neuron
    table.insert(network.neurons,tempN)
		network.neurons[i] = tempN
    table.insert(network.weightIndex,weightTracker)
    weightTracker = weightTracker + 1
   -- print("neuron input value: "..network.neurons[i].value)
	end

	--create outputs
	for o=1,#outputs do
    tempO = newNeuron()
    tempO.value =  math.random(500,800)
    tempO.weightIndex = network.weightIndex
    table.insert(network.neurons,tempO)
		network.neurons[MaxNodes+o] = tempO
    -- print("neuron output value: "..network.neurons[o].value)
	end


  print("network size after adding inputs and outputs value"..#network.neurons)
  --add neurons to genome network
  for i = 1, #network.neurons do
    table.insert(genome.network,network.neurons[i])
  end

  --create connection gene to match inputs index just for the output
  for i = 1, #inputs do
    connN = connectionGene()
    connN.input = genome.network[i].value
    connN.out = network.neurons[MaxNodes+1]
    table.insert(network.weights,connN)
    end

  --copy the connection genes
  for i = 1, #network.weights do
    table.insert(genome.genes,network.weights[i])
    end

  return genome
end


function evaluateNetwork(genome)
  --obtain all connections to a node and shit out output
  for i = 1, #genome.network do
    tempW = {} --all associated
    tempWout = {} --just the outs
    oldValue = genome.network[i].value
    --update node AND gene
    --obtain node value
    --TempNode = genome.network[i]
    print("genome network neuron " ..genome.network[i].value)
    print("weight index values " ..#genome.network[i].weightIndex)
    --obtain weights through index
    for j = 1, #genome.network[i].weightIndex do
      --grab index and loop through network weights to check if its there
      tempWeightIndx = genome.network[i].weightIndex[j]
      print("weight index to look for"..tempWeightIndx)
      for k = 1, #genome.genes do
        if tempWeightIndx == k then
          --clear a table: for k,v in pairs(tab) do tab[k]=nil end
          --add to tempW
          print("Associated gene in"..genome.genes[k].input)
          print("Associated gene in"..genome.genes[k].out.value)
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
        if tempW[b].out.value == genome.network[i].value and tempW[b].status == true then
          table.insert(tempWout, tempW[b])
        end
    end
    print("found refined connection genes : "..#tempWout)
    --loop through all cases of tempWout, obtain ins, multiply by weights get new out value and REPLACE the gene with the new tempWouts
    sum = 0
    activation = 0
    for m = 1, #tempWout do
      print("input value from filtered "..tempWout[m].input)
      sum = sum + (tempWout[m].input * tempWout[m].weight)
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
          if genome.genes[q].innovation == tempWout[p].innovation and genome.genes[q].out.value == oldValue then
            print("previous "..genome.genes[q].out.value)
            genome.genes[q]= tempWout[p]
            print("current "..genome.genes[q].out.value)
          end
          --if same innovation number and same input
          if genome.genes[q].innovation == tempWout[p].innovation and genome.genes[q].input == oldValue then
            print("previous "..genome.genes[q].input)
            genome.genes[q]= tempWout[p]
            print("current "..genome.genes[q].input)
          end

          end
        end
      --REPLACE NODE
    end
    print("new genome value".. genome.network[i].value)
  end
  end
 testPropagate = createNewGenome()
 testPropagate = createNewNetwork(testPropagate)
 f_to_pay_respects = evaluateNetwork(testPropagate)
 --print("gene value"..testPropagate.network[2].value)
 --print("gene out value"..testPropagate.genes[2].out.value)
 --I want to find out how many weight values are attatched to me




--  function newNeuron()
--	local neuron = {}
--  neuron.weightIndex = {} --keeps track of attatched weights (matches with node array index
--  neuron.inStatus = {} --0 if an input 1 if an output (matches with node array index
--	neuron.incoming = {} --data from previous thingis
--	neuron.value = 0.0 -- current neuron value
--	return neuron
--end
