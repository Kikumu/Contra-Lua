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
--Now for de crossover of my little croissants

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



function mutateConnectionGene(gene)
--take 2 random nodes and adds a connection between them if none are available

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
if state1 == 1 then
connectionGeneMutate1 = connectionGene()
connectionGeneMutate1.input = node1.input
connectionGeneMutate1.innovation = node1.innovation + 1
connectionGeneMutate1.out = node2.input
--connectionGeneMutate1.output = a(i*w)?
connectionGeneMutate1.out= node2.input
table.insert(gene,connectionGeneMutate1)
end

--if state2 == 1 then

--end



print('state1: '..state1)
print('state2: '..state2)
return gene
end


function mutateNodeGene(gene)
--pick a random node/neuron/connection gene
genepos = math.random(1,#gene)
node1 = gene[genepos]


--disable gene
node1.status = false

--put it back

gene[genepos] = node1

print("node1 in: "..node1.input)
--node one in to new connection to new neuron

connectionGeneMutate1 = connectionGene()

connectionGeneMutate1.input = node1.input
connectionGeneMutate1.innovation = node1.innovation + 1
--connectionGeneMutate1.output = a(i*w)?
--add neuron to a table
incrementneurons = #neurons + 1
neurons[incrementneurons] = math.random()
connectionGeneMutate1.output = neurons[incrementneurons]
table.insert(gene,connectionGeneMutate1)

connectionGeneMutate2 = connectionGene()
connectionGeneMutate2.input = connectionGeneMutate1.out
connectionGeneMutate2.innovation = connectionGeneMutate1.innovation + 1

--connectionGeneMutate2.output = a(i*w)?
connectionGeneMutate2.out = node1.out
table.insert(gene,connectionGeneMutate2)

return gene
end



--for crossover purposes(Goal is to add to g1)
function matchingGenes(g1,g2)
matchedGenes ={}
g1Length = #g1
g2Length = #g2

for i = 0, g1Length do
for j = 0, g2Length do
if g1[i].innovation == g2[j].innovation then
if selectionChance > math.random()
table.insert(matchedGenes,g1[i])
else
table.insert(matchedGenes,g2[j])
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
g2MaxInnovation = g2[#g2].innovation
g1MaxInnovation = g1[#g1].innovation
g1Length = #g1
g2Length = #g2

for i = 0, g1Length do
for j = 0, g2Length do
if g1[i].innovation ~= g2[j].innovation and g2[j].innovation < g1MaxInnovation then
table.insert(disjointedGenes, g2[i])
end
end
end


return disjointedGenes

end


--confident(not tested)
function excessGenes(g1,g2)
ExcessGenesTable = {}
--find excess genes in g2(which are located in g1)
maxInnovationg1 = g1[#g1].innovation
maxInnovationg2 = g2[#g2].innovation
g1Length = #g1
g2Length = #g2

if maxInnovationg1 > maxInnovationg2
for i = 0, g1Length do
if g1[i].innovation > maxInnovationg2 then
table.insert(ExcessGenesTable, g1[i])
end
end
end

if maxInnovationg2 > maxInnovationg1
for i = 0, g2Length do
if g2[i].innovation > maxInnovationg1 then
table.insert(ExcessGenesTable, g2[i])
end
end
end
return ExcessGenesTable
end


--sorted child gene according to innovation number
function selectionSort(DisjointedGenesArr, ExcessGenesArr, MatchingGenesArr)

genesToSort = {}
table.insert(genesToSort,DisjointedGenesArr)
table.insert(genesToSort,ExcessGenesArr)
table.insert(genesToSort,MatchingGenesArr)

for i = 0, #genesToSort do
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

return genesToSort
end


--breed
function crossover(genome) --genome in the sense that you are passing to this function a set of genomes

childGenome= {}

childlength = 0
gene1 = genome[math.random(1,#genome)]
gene2 = genome[math.random(1,#genome)]



end

--if they are not equal, it could either be: an excess gene OR a disjoint gene.
--how do i tell the difference....
--an excess gene exceeds parameters(innovation number) of another gene
if gene1[i].innovation ~= gene2[i].innovation then

end
--[[what about excesses and disjoints. Remember, Excesses are above the innovation numbers
of a genome while disjoints are innovation numbers in which one parent has and the other
hasnt in the middle]]

--lets start with disjoints
--given 2 genomes





end

end






function species(genes)

speciesPool = {}

return speciesPool
end


x = newGenome()
y = mutateNodeGene(x)


print("number is "..#inputs)
print(x[2].weight)
print(#neurons)
print(y[5].innovation)
print(neurons[2])




