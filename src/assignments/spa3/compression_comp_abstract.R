svd = function(x){
  512 * x + x * x + x * 512
} 

koh = function(nodes, block.size){
    nodes*block.size + 512**2 / block.size
}

nodes      = seq(1,100, length.out=100)
block.size = seq(1,100, length.out=100)
z = outer(nodes,block.size, "*")

y = 512**2/block.size
d = outer(z,y, "+")
persp(nodes,block.size,d, col='blue')
# 

