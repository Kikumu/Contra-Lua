mutationChance = 0.5
geneActivationChance = 0.5

--species classification variables

disjointment = 0.2 --c1
excessGenes = 0.3 --c2
weightImportance = 0.1 --c3

--E is number of disjointed connection genes
--D is number of excess connection genes
--W is weight value
--N number of connection genes

--general formula on paper: ( c1*E/N + c2*D/N + c3*W)


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
function newGene()

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

node1 = gene[math.random(1,3)]
node2 = gene[math.random(1,3)]
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



print('state1: '..state1)
print('state2: '..state2)
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
connectionGeneMutate1.output = 1
table.insert(gene,connectionGeneMutate1)

connectionGeneMutate2 = connectionGene()
connectionGeneMutate2.inputs = connectionGeneMutate1.out
connectionGeneMutate2.innovation = connectionGeneMutate1.innovation + 1

--connectionGeneMutate2.output = a(i*w)?
connectionGeneMutate2.output = node1.out
table.insert(gene,connectionGeneMutate1)

return gene
end

function crossover(gene1,gene2)

end

function species(genes)

speciesPool = {}


return speciesPool
end


x = newGene()
y = mutateNodeGene(x)

print("number is "..#inputs)
print(x[2].weight)
print(#y)



