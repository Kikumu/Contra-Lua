mutationChance = 0.5
geneActivationChance = 0.5
selectionChance = 0.3 --this value is used to choose between 2 connection genes of the same innovation
--species classification variables

disjointment = 0.2 --c1
excessGenes = 0.3 --c2
weightImportance = 0.1 --c3

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
1) Take it out for a test spin to test out functions created
2) Find a way to keep track of the fitness of each genome
3) Create function to organise a group of genes into different speciesPool

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
    return 1.0 / (1.0 + exp(-x))
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

