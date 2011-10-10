from django.db import models
from math import floor

# Create your models here.
#primer TEXT, box TEXT, row TEXT, col INT, seq TEXT,  date TEXT , comments TEXT)')

class Order(models.Model):
    norder=models.CharField(max_length=30)
    date=models.DateField()
    
    def __unicode__(self):
        return u'%s %s' % (self.norder, self.date)
    
class Pending(models.Model):
    
    name=models.CharField(max_length=30)
    seq=models.CharField(max_length=100)
    #order= models.ForeignKey(Order)
    def __unicode__(self):
        return self.name

class Primer(models.Model):
    
    name = models.CharField(max_length=30)
    box = models.IntegerField(max_length=2)
    row = models.CharField(max_length=1)
    col = models.IntegerField(max_length=2)
    seq = models.CharField(max_length=100)
    #order= models.ForeignKey(Order)
    comment = models.CharField(max_length=100)
    loc_id = models.IntegerField(max_length=10)
    
    def __unicode__(self):
        return self.name
    
    @staticmethod
    def get_loc_id(box, row, col):
        return (box-1)*81 + (ord(row)-ord('A'))*9 + (col-1)
    
    @staticmethod
    def loc_id2location(loc_id):
        box = int(floor(loc_id / 81)) + 1
        residual = loc_id % 81
        row = chr(int(floor(residual / 9)) + 65)
        residual = residual % 9
        col = residual + 1
        return box, row, col
    