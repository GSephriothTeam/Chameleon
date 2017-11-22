# Chameleon
CS550 &amp; CS554 course project

服务器配置
- node1 23.101.168.192
- node2 52.162.241.179
- node3 52.237.158.174
  
三个节点配置相同，将以下配置中的1，2，3相互替换即可
- root用户 用户名:node1 密码:Sparkhadoopnode1
- Hadoop和Spark用户 用户么:hduser 密码:123456

####!!Hadoop version: 2.6.5

####!!Spark version: 2.0.2


代码源文件位于 "src/" 文件夹下
可执行文件将被编译到 “bin／” 文件夹下
默认使用Ant编译，编译配置文件为 “build.xml”

“src/Chameleon” 为Chameleon算法代码，包含csv文件读取，运行结果将被保存到res.txt
"src/Data" 为测试数据生成代码
"bin/Data1.txt"和"src/graphData.txt" 用于测试程序能否正常运行
"bin/Data*.txt" 用于测试程序的性能

####修改编译时的主函数:
58行 <mainClass>Spark.WordCount</mainClass> 指定主函数为 WordCount
59行 <mainClass>Chameleon.Client</mainClass> 指定主函数为 Client