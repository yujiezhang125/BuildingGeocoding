# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 10:22:45 2020

@author: yujiezhang125
"""


import numpy as np
import pandas as pd
import geopandas as gpd
import math
import requests
from requests.exceptions import ReadTimeout, ConnectTimeout
import json
import time
pd.set_option('display.max_columns', 500)

# converting coords
class LngLatTransfer():

    def __init__(self):
        self.x_pi = 3.14159265358979324 * 3000.0 / 180.0
        self.pi = math.pi  # π
        self.a = 6378245.0  # 长半轴
        self.es = 0.00669342162296594323  # 偏心率平方
        pass

    def GCJ02_to_BD09(self, gcj_lng, gcj_lat):
        """
        实现GCJ02向BD09坐标系的转换
        :param lng: GCJ02坐标系下的经度
        :param lat: GCJ02坐标系下的纬度
        :return: 转换后的BD09下经纬度
        """
        z = math.sqrt(gcj_lng * gcj_lng + gcj_lat * gcj_lat) + 0.00002 * math.sin(gcj_lat * self.x_pi)
        theta = math.atan2(gcj_lat, gcj_lng) + 0.000003 * math.cos(gcj_lng * self.x_pi)
        bd_lng = z * math.cos(theta) + 0.0065
        bd_lat = z * math.sin(theta) + 0.006
        return bd_lng, bd_lat

    def BD09_to_GCJ02(self, bd_lng, bd_lat):
        '''
        实现BD09坐标系向GCJ02坐标系的转换
        :param bd_lng: BD09坐标系下的经度
        :param bd_lat: BD09坐标系下的纬度
        :return: 转换后的GCJ02下经纬度
        '''
        x = bd_lng - 0.0065
        y = bd_lat - 0.006
        z = math.sqrt(x * x + y * y) - 0.00002 * math.sin(y * self.x_pi)
        theta = math.atan2(y, x) - 0.000003 * math.cos(x * self.x_pi)
        gcj_lng = z * math.cos(theta)
        gcj_lat = z * math.sin(theta)
        return gcj_lng, gcj_lat

    def WGS84_to_GCJ02(self, lng, lat):
        '''
        实现WGS84坐标系向GCJ02坐标系的转换
        :param lng: WGS84坐标系下的经度
        :param lat: WGS84坐标系下的纬度
        :return: 转换后的GCJ02下经纬度
        '''
        dlat = self._transformlat(lng - 105.0, lat - 35.0)
        dlng = self._transformlng(lng - 105.0, lat - 35.0)
        radlat = lat / 180.0 * self.pi
        magic = math.sin(radlat)
        magic = 1 - self.es * magic * magic
        sqrtmagic = math.sqrt(magic)
        dlat = (dlat * 180.0) / ((self.a * (1 - self.es)) / (magic * sqrtmagic) * self.pi)
        dlng = (dlng * 180.0) / (self.a / sqrtmagic * math.cos(radlat) * self.pi)
        gcj_lng = lat + dlat
        gcj_lat = lng + dlng
        return gcj_lng, gcj_lat

    def GCJ02_to_WGS84(self, gcj_lng, gcj_lat):
        '''
        实现GCJ02坐标系向WGS84坐标系的转换
        :param gcj_lng: GCJ02坐标系下的经度
        :param gcj_lat: GCJ02坐标系下的纬度
        :return: 转换后的WGS84下经纬度
        '''
        dlat = self._transformlat(gcj_lng - 105.0, gcj_lat - 35.0)
        dlng = self._transformlng(gcj_lng - 105.0, gcj_lat - 35.0)
        radlat = gcj_lat / 180.0 * self.pi
        magic = math.sin(radlat)
        magic = 1 - self.es * magic * magic
        sqrtmagic = math.sqrt(magic)
        dlat = (dlat * 180.0) / ((self.a * (1 - self.es)) / (magic * sqrtmagic) * self.pi)
        dlng = (dlng * 180.0) / (self.a / sqrtmagic * math.cos(radlat) * self.pi)
        mglat = gcj_lat + dlat
        mglng = gcj_lng + dlng
        lng = gcj_lng * 2 - mglng
        lat = gcj_lat * 2 - mglat
        return lng, lat

    def BD09_to_WGS84(self, bd_lng, bd_lat):
        '''
        实现BD09坐标系向WGS84坐标系的转换
        :param bd_lng: BD09坐标系下的经度
        :param bd_lat: BD09坐标系下的纬度
        :return: 转换后的WGS84下经纬度
        '''
        lng, lat = self.BD09_to_GCJ02(bd_lng, bd_lat)
        return self.GCJ02_to_WGS84(lng, lat)

    def WGS84_to_BD09(self, lng, lat):
        '''
        实现WGS84坐标系向BD09坐标系的转换
        :param lng: WGS84坐标系下的经度
        :param lat: WGS84坐标系下的纬度
        :return: 转换后的BD09下经纬度
        '''
        lng, lat = self.WGS84_to_GCJ02(lng, lat)
        return self.GCJ02_to_BD09(lng, lat)

    def _transformlat(self, lng, lat):
        ret = -100.0 + 2.0 * lng + 3.0 * lat + 0.2 * lat * lat + 0.1 * lng * lat + 0.2 * math.sqrt(math.fabs(lng))
        ret += (20.0 * math.sin(6.0 * lng * self.pi) + 20.0 *
                math.sin(2.0 * lng * self.pi)) * 2.0 / 3.0
        ret += (20.0 * math.sin(lat * self.pi) + 40.0 *
                math.sin(lat / 3.0 * self.pi)) * 2.0 / 3.0
        ret += (160.0 * math.sin(lat / 12.0 * self.pi) + 320 *
                math.sin(lat * self.pi / 30.0)) * 2.0 / 3.0
        return ret

    def _transformlng(self, lng, lat):
        ret = 300.0 + lng + 2.0 * lat + 0.1 * lng * lng + 0.1 * lng * lat + 0.1 * math.sqrt(math.fabs(lng))
        ret += (20.0 * math.sin(6.0 * lng * self.pi) + 20.0 *
                math.sin(2.0 * lng * self.pi)) * 2.0 / 3.0
        ret += (20.0 * math.sin(lng * self.pi) + 40.0 *
                math.sin(lng / 3.0 * self.pi)) * 2.0 / 3.0
        ret += (150.0 * math.sin(lng / 12.0 * self.pi) + 300.0 *
                math.sin(lng / 30.0 * self.pi)) * 2.0 / 3.0
        return ret

    def WGS84_to_WebMercator(self, lng, lat):
        '''
        实现WGS84向web墨卡托的转换
        :param lng: WGS84经度
        :param lat: WGS84纬度
        :return: 转换后的web墨卡托坐标
        '''
        x = lng * 20037508.342789 / 180
        y = math.log(math.tan((90 + lat) * self.pi / 360)) / (self.pi / 180)
        y = y * 20037508.34789 / 180
        return x, y

    def WebMercator_to_WGS84(self, x, y):
        '''
        实现web墨卡托向WGS84的转换
        :param x: web墨卡托x坐标
        :param y: web墨卡托y坐标
        :return: 转换后的WGS84经纬度
        '''
        lng = x / 20037508.34 * 180
        lat = y / 20037508.34 * 180
        lat = 180 / self.pi * (2 * math.atan(math.exp(lat * self.pi / 180)) - self.pi / 2)
        return lng, lat

# # Gaode API
def gaode(i,data,amap_dataframe):
    para = {
        'address': "重庆" + str(data.iloc[i,:].位置) + str(data.iloc[i,:].传统商贸批发市场名称), ##################
        'key':'4c65d8ba694fe67466d28fe2e29d50ee'
    }
    url = 'https://restapi.amap.com/v3/geocode/geo?'
    
    #print(r.json())
    print(i)
    try:
        r = requests.get(url, para, timeout=30)
        if r.json()['info']=='OK':
            if len(r.json()['geocodes'])==1:
                print(r.json()['geocodes'][0]['level'])

                amap_dataframe.at[i,'geocode'] = r.json()['geocodes'][0]['location']
                amap_dataframe.at[i,'formatted_address'] = r.json()['geocodes'][0]['formatted_address']
                amap_dataframe.at[i,'ID'] = data.iloc[i,:].ID
                amap_dataframe.at[i,'level'] = r.json()['geocodes'][0]['level']
                # amap_dataframe.at[i,'address'] = data.iloc[i,:].address
            
            else:
                amap_dataframe.at[i,'geocode'] = 'null'
                amap_dataframe.at[i,'formatted_address'] = 'null'
                amap_dataframe.at[i,'ID'] = data.iloc[i,:].ID
                amap_dataframe.at[i,'level'] = 'null'
                # amap_dataframe.at[i,'address'] = data.iloc[i,:].address
    except (ReadTimeout, ConnectTimeout):
        pass

# loop geocoding function
def loopgeo(data, col_names):
    # create empty df to save api results
    amap0 = pd.DataFrame(columns=col_names)

    for i in range(len(data)):
        gaode(i, data, amap0)

    amap0['lon']= amap0.geocode.str.split(",", expand=True)[0]
    amap0['lat']= amap0.geocode.str.split(",", expand=True)[1]

    amap0['lon84']=0.0
    amap0['lat84']=0.0
    
    t = LngLatTransfer()
    
    for i, row in amap0.iterrows():
        if row.geocode!='null':
            gcj_lat = float(row.lat)
            gcj_lng = float(row.lon)
            amap0.loc[i,'lon84'], amap0.loc[i,'lat84'] = t.GCJ02_to_WGS84(gcj_lng, gcj_lat)

    result0 = data.merge(amap0, on = 'ID')
    return result0


# # baidu API
def baidu(i,data,amap_dataframe):
    para = {
        'output':'json',
        'address': "重庆" + str(data.iloc[i,:].位置) + str(data.iloc[i,:].传统商贸批发市场名称), ####################
        'ak':'rAGUuvw74LvPKuNiwsTkAMCRRVeuwppk'
    }
    url = 'http://api.map.baidu.com/geocoding/v3/?'
    print(i)
    
    try:
        r = requests.get(url, para, timeout=30)
        print(r.json())
        if r.json()['status']==0:
            if (len(r.json()['result']['location'])==2): #& (r.json()['result']['level'] not in ['省','城市','区县','乡镇','村庄','道路'])       

                amap_dataframe.at[i,'lat'] = r.json()['result']['location']['lat']
                amap_dataframe.at[i,'lon'] = r.json()['result']['location']['lng']
                amap_dataframe.at[i,'level'] = r.json()['result']['level']

    except (ReadTimeout, ConnectTimeout):
        amap_dataframe.at[i,'geocode'] = 'null'
        pass

def loop_baidu(data): 
    for i in range(len(data)):
        baidu(i, data, data)
    
    t = LngLatTransfer()
    
    data['lon84']=0.0
    data['lat84']=0.0
    for i, row in data.iterrows():
        gcj_lat = float(row.lat)
        gcj_lng = float(row.lon)
        data.loc[i,'lon84'], data.loc[i,'lat84'] = t.BD09_to_WGS84(gcj_lng, gcj_lat)
        
        
# # tencent API
def tencent(i,data,amap_dataframe):
    para = {
        'address': "重庆" + str(data.iloc[i,:].位置) + str(data.iloc[i,:].传统商贸批发市场名称), ###################
        'key':'774BZ-BPUKX-4IP4R-ZYLX4-7W5C2-J4BSH'
    }
    url = 'https://apis.map.qq.com/ws/geocoder/v1/?'
    print(i)
    
    try:
        r = requests.get(url, para, timeout=30)
        print(r.json())
        if r.json()['status']==0:
            if (len(r.json()['result']['location'])==2): #& (r.json()['result']['level'] not in ['省','城市','区县','乡镇','村庄','道路'])       

                amap_dataframe.at[i,'lat'] = r.json()['result']['location']['lat']
                amap_dataframe.at[i,'lon'] = r.json()['result']['location']['lng']
                amap_dataframe.at[i,'level'] = r.json()['result']['level']
                amap_dataframe.at[i,'formatted_address'] = r.json()['result']['title']
                time.sleep(0.2)
            
    except (ReadTimeout, ConnectTimeout):
        amap_dataframe.at[i,'geocode'] = 'null'
        pass

def loop_tencent(data):
    for i in range(len(data)):
        tencent(i, data, data)
    
    t = LngLatTransfer()
    
    data['lon84']=0.0
    data['lat84']=0.0
    for i, row in data.iterrows():
        gcj_lat = float(row.lat)
        gcj_lng = float(row.lon)
        data.loc[i,'lon84'], data.loc[i,'lat84'] = t.GCJ02_to_WGS84(gcj_lng, gcj_lat)

########### import data
mkt = pd.read_excel("现状传统商贸批发市场名录和会展中心名录.xlsx")
mkt = mkt.rename(columns={'序号':'ID'})

################# geocoding!
df = mkt
filter1 = ['null','道路','省', '市', '区县', '未知', '乡镇']#,'村庄'
filter2 = ['区县','NoClass','城市']#'道路','村庄','乡镇','',
filter3 = [1,2,7]

print("geocoding AMAP...")
geocode_gaode = loopgeo(df, ['ID','geocode','formatted_address', 'level'])
gaoder = geocode_gaode[~geocode_gaode.level.isin(filter1)].reset_index(drop=True)

if len(geocode_gaode[geocode_gaode.level.isin(filter1)]) > 0:  
    geocode_baidu = geocode_gaode[geocode_gaode.level.isin(filter1)].reset_index(drop=True)
    print("geocoding BAIDU...")
    loop_baidu(geocode_baidu)
    baidur = geocode_baidu[~geocode_baidu.level.isin(filter2)].reset_index(drop=True)
    
    if len(geocode_baidu[geocode_baidu.level.isin(filter2)]) > 0:     
        geocode_qq = geocode_baidu[geocode_baidu.level.isin(filter2)].reset_index(drop=True)
        print("geocoding TENCENT...")
        loop_tencent(geocode_qq)
        qq = geocode_qq[~geocode_qq.level.isin(filter3)].reset_index(drop=True)

geocode_result = pd.concat([pd.concat([qq, baidur]),gaoder]).reset_index(drop=True)

############# export
gaoder.to_csv("geocode_result_exhibit.csv",encoding='gb18030')
geocode_result.to_csv("geocode_result_sport2.csv",encoding='gb18030')


