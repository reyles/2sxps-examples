#Downloads UVOT data for a specific position.
from __future__ import division
import numpy as np
import time
import warnings
from numpy import *
from bs4 import BeautifulSoup
import requests
import selenium
# Using Chrome to access web
from selenium import webdriver
from selenium.webdriver import Chrome
import os

#Enter source name and it makes up the search string
src = 'J105534.3-012614'

#Makes the directory for depositing the data
os.system('mkdir ./data/%s' % src)

#Sets up webdriver
options = webdriver.ChromeOptions()
options.add_argument('--ignore-certificate-errors')
options.add_argument('--incognito')
#This line commented out so the download can be monitored
#options.add_argument('--headless')
options.add_argument('--disable-extensions')
options.add_experimental_option("detach", True)
prefs = {'download.default_directory' : '/Users/raje1/Documents/Projects/2SXPS/uvot/data/%s' % src}
options.add_experimental_option('prefs', prefs)
global driver
driver = Chrome(options = options)

#Cone search - Swift webpage has standardised address
driver.get('https://www.swift.ac.uk/swift_live/doSimpleSearch.php?catname=swiftmastr&searchpos=%s%%3A%s%%3A%s+%s%%3A%s%%3A%s&searchrad=12&searchepoch=2000&getWhat=1&sortBy=_r&dandr=on&coordType=0&retEpoch=2000&retType=0&retWhere=1'
           % (src[1:3],src[3:5],src[5:9],src[9:12],src[12:14],src[14:16]))
#Required xpath for objects
xpath_obj = '/html/body/div/div[3]/table/tbody/tr[1]/td[1]/a[2]'
driver.find_element_by_xpath(xpath_obj).click()
#Required xpath for download
xpath_obj = '/html/body/div/table/tbody/tr[last()]/td[2]/a'
driver.find_element_by_xpath(xpath_obj).click()
#Required xpaths to only download uvot data
xpath_auxil = '/html/body/div/div[4]/form/div[1]/table/tbody/tr[1]/td[1]/input'
xpath_uvot = '/html/body/div/div[4]/form/div[1]/table/tbody/tr[4]/td[1]/input'
xpath_download = '/html/body/div/div[4]/form/div[1]/input[3]'
driver.find_element_by_xpath(xpath_auxil).click()
driver.find_element_by_xpath(xpath_uvot).click()
driver.find_element_by_xpath(xpath_download).click()
#While loop and sleep to wait for download to be completed
files = os.listdir('data/%s' % src)
while 'download.tar' not in files:
    files = os.listdir('data/%s' % src)
time.sleep(5)
#Extracts files
os.system('tar -xvf data/%s/download.tar -C /Users/raje1/Documents/Projects/2SXPS/uvot/data/%s' % (src,src))
os.system('rm data/%s/download.tar' % src)

print('Done!')
