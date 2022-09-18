import datetime
import os
import random
import re
import shutil
import time
import zipfile

import numpy as np
import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By

current_path = os.getcwd()

# Start time
t1 = datetime.datetime.now()
print('Start time:', t1.strftime('%Y-%m-%d %H:%M:%S'))

chrome_options = Options()
chrome_options.add_argument('--headless')
chrome_options.add_argument('--disable-gpu')
headers = 'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/96.0.4664.9 Safari/537.36'

# # test
# target  = 'VEGFR'   # target to download
# seconds = 20        # rest time


import PySimpleGUI as sg

theme_name_list = sg.theme_list()  # theme list
# sg.theme_previewer()            # to show theme
sg.theme(random.choice(theme_name_list))  # set theme
layout = [
    [sg.Text('Please Enter the target name           '), sg.InputText()],
    [sg.Text('Please Enter the sleep time seconds'), sg.InputText()],
    [sg.Button('Ok'), sg.Button('Cancel')]
]
# create window
window = sg.Window('download inhibitors of any target', layout, font=("Arial", 20), default_element_size=(50, 1),
                   size=(1000, 500))
# Event loop and get values
while True:
    event, values = window.read()
    if event in (None, 'Ok'):
        break
    print('You entered ', values[0])

target = values[0]  # target to download
seconds = int(values[1])  # rest time
window.close()

try:
    os.mkdir(target)
except FileExistsError:
    pass

driver_list = [i for i in os.listdir() if re.search(r"driver.exe", i)]
driver_list.reverse()
driver_name = driver_list[0]

if driver_name == "msedgedriver.exe":
    driver = webdriver.Edge(service=Service(r"msedgedriver.exe"))
elif driver_name == "chromedriver.exe":
    driver = webdriver.Chrome(service=Service(r"chromedriver.exe"))
else:
    print('you should add a driver in the current folder according to your version of broswer')

url = 'https://www.ebi.ac.uk/chembl/g/#browse/activities'
driver.get(url)
driver.maximize_window()  # 页面最大化
time.sleep(seconds)
# 获取该输入框的ID  //*[@id="ESActivitity-Table"]/div/div[2]/div[1]/div[3]/div/div/input
input_element = driver.find_element(By.XPATH, '//*[@id="ESActivitity-Table"]/div/div[2]/div[1]/div[3]/div/div/input')
# 清除该输入框中的原本内容
input_element.clear()
# 向该输入框中添加搜索词
input_element.send_keys("{0}".format(target))
# 模拟点击
button_element = driver.find_element(By.XPATH,
                                     '//*[@id="ESActivitity-Table"]/div/div[2]/div[1]/div[3]/div/div/div[1]/a')

# 点击展开csv
button_element.click()
time.sleep(seconds)

# 点击csv展开下载选项
csv_btn = driver.find_element(By.XPATH, '//*[@id="GladosMainContent"]/div/div/div[1]/div[2]/div[1]/div[3]/div/a[1]')
csv_btn.click()

# Generating Download File (100%)
time.sleep(seconds)

click_here = driver.find_element(By.XPATH,
                                 '//*[@id="GladosMainContent"]/div/div/div[1]/div[2]/div[2]/div/div/div[4]/a[1]')
click_here.click()  # ZIP文件 默认下载到 C:\Users\Administrator\Downloads

file_name = click_here.get_attribute('href').split('/')[-1:][0]

time.sleep(seconds)

# default absolute download path
import getpass

root_user = getpass.getuser()
download_path = "C:\\Users\\{}\\Downloads".format(root_user)
srcfile = os.path.join(download_path, file_name)

# file_list2 = []
# for root, dirs, files in os.walk(download_path):
#     sorted(files, key=lambda x: os.path.getmtime(os.path.join(dirs, x)))
#     for file in files:
#         matchObj = re.match(r'~$', file, re.I)
#         if not matchObj:
#             if os.path.splitext(file)[1] in [".zip"]:
#                 file_list2.append(os.path.join(root, file))
#
# srcfile = file_list2[0]
# copy file
newfile = srcfile.replace(download_path, current_path + '\\' + target)
shutil.copyfile(srcfile, newfile)
time.sleep(seconds)
zipObj = zipfile.ZipFile(newfile)
zip_list = zipObj.namelist()
# # 循环解压文件到指定绝对路径
for f in zip_list:
    zipObj.extract(f, current_path + '\\' + target)
zipObj.close()  # 关闭文件，必须有，释放内存

# 第一步：导入数据
file = zip_list[0]
df = pd.read_csv(target + '/' + file, sep=';')
# 筛选数据，保留感兴趣的记录
df1 = df.loc[df["Standard Relation"] == "'='", :]
df1 = df1.loc[df1["Standard Value"] >= 0, :]
# 筛选数据，保留感兴趣的变量
df2 = df1.loc[:, ['Molecule ChEMBL ID', 'Smiles', 'Standard Value']]

f = open(target + '/{0}.txt'.format(target), mode='w', encoding='utf-8')
smi_list = df2['Smiles'].tolist()
for smi in smi_list:
    print(smi)
    f.write(str(smi) + '\n')
f.close()


# 定义计算函数，衍生新的变量 如: activity pIC50
def classifier(x):
    if x > 1000:
        return 0
    elif 100 <= x <= 1000:
        return 1
    elif x < 100:
        return 2


def calcpIC50(x):
    return -np.log10(x * 10 ** (-9))


df2['activity'] = df2['Standard Value'].apply(lambda x: classifier(x))
df2['pIC50'] = df2['Standard Value'].apply(lambda x: calcpIC50(x))
df3 = df2.sort_values(by=["Standard Value"], ascending=True)
df3.index = range(len(df3))
df3.to_excel(target + '/{}.xlsx'.format(target), index=False)

print('Start time:', t1.strftime('%Y-%m-%d %H:%M:%S'))
# End time
t2 = datetime.datetime.now()
print('End time:', t2.strftime('%Y-%m-%d %H:%M:%S'))
delta = t2 - t1
if delta.seconds > 3600:
    print('Time Consuming: '
          + str(delta.seconds // 3600) + ' h '
          + str(delta.seconds % 3600 // 60) + ' min '
          + str(delta.seconds % 60) + ' seconds ')
elif delta.seconds > 60:
    print('Time Consuming: '
          + str(delta.seconds // 60) + ' min '
          + str(delta.seconds % 60) + ' seconds ')
else:
    print('Time Consuming: ' + str(delta.seconds) + ' seconds')

print("finished!")
