/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CBCTop;

import jsat.classifiers.DataPoint;

/**
 *
 * @author dell
 */
public class PairResiduIdFeature {
    private String residuId;
    private String residuName;
    private DataPoint feature;
   // private 

    public PairResiduIdFeature(String id,String name, DataPoint p){
    this.residuId=id;
    this.residuName = name;
    this.feature = p;
    }

    public DataPoint getFeature() {
        return feature;
    }

    public String getResiduId() {
        return residuId;
    }
    
}
