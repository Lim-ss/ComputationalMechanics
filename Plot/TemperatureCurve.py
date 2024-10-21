import matplotlib.pyplot as plt

# 读取文件中的数据
data_file = './../ComputationalMechanics/ComputationalMechanics/Temperature.data'  # 文件路径
x_values = []  # 存储 x 值
y_values = []  # 存储 y 值

with open(data_file, 'r') as file:
    for line in file:
        x, y = line.split()  # 假设文件中的数据以空格分隔
        x_values.append(float(x))
        y_values.append(float(y))

# 绘制曲线图
plt.plot(x_values, y_values)
plt.xlabel('x(m)')
plt.ylabel('T(K)')
plt.title('Temperature curve')
plt.show()