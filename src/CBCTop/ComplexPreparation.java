/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package CBCTop;

import ComplexFeature.AAIndex;
import ComplexFeature.PSAI;
import ComplexStructure.Aminoacid;
import ComplexStructure.Complex;
import ExtractDataTraining.ExtractLogOdds;
import ExtractDataTraining.ExtractProtusionIndex;
import ExtractDataTraining.ExtractSphereExposure;
import ExtractFromFile.IO;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import jsat.classifiers.DataPoint;
import jsat.linear.DenseVector;

/**
 *
 * @author dell
 */
public class ComplexPreparation {
    private Complex complex;
    private String pdbid;
    private String chain;
    private double treshold;
    public ComplexPreparation(String pdbid, String c, double t){
    complex = new Complex(pdbid,c);
    this.pdbid = pdbid;
    this.chain = c;
    this.treshold=t;
    
    }
    public ComplexPreparation(String path,String pdbid, String c, double t){
        complex = new Complex(path,pdbid,c);
        this.pdbid = pdbid;
        this.chain = c;
        this.treshold=t;
    }
    public ArrayList<PairResiduIdFeature> getPairIdNFeature(String path,String akhiran) throws IOException{
    //bangun struktur complex dari file pdb
    // set parameter psaia
   // this.setExposedResidueFromPSAIA();
   this.setExposedRFromPSAIA(pdbid, path, akhiran);
    //ambil id residu terekspose 
    //ekstraksi ciri residu terekspose
    this.setLogOdd();
    this.setSphereExposure();
    ArrayList<Aminoacid> listAA = complex.getAminoacidVector();
    ArrayList<PairResiduIdFeature> idNF= new ArrayList();
    for(Aminoacid aa:listAA){
        
        double  rsa = aa.getRsaTien2013();
        if(rsa>this.treshold){
            List<Double> atrib = new ArrayList();
            atrib.add(aa.getRsaTien2013());
            
            atrib.addAll(aa.getpASA().getParamASAD());
            atrib.add(aa.getCaBfactor());
            atrib.add(aa.getBfactor());
            atrib.add(aa.getLo());
            atrib.addAll(aa.getPsai().getListPSAI());
            atrib.addAll(this.getAAIList(aa));
            DenseVector vec = new DenseVector(atrib);
            DataPoint p= new DataPoint(vec);
            PairResiduIdFeature pridf = new PairResiduIdFeature(aa.getresidueID(),aa.getresidueName(),p);
            idNF.add(pridf);
            
        }
        
    }
    System.out.println("idNF"+ idNF.size());
    return idNF;
    //lakukan klasifikasi level residu
    
    //System.out.println("jumlah pasangan residu dan feature"+ idNF.size());
    //temukan klaster-klaster berdasarkan kedekatan spasial
    //temukan kandidat terbaik
    }

    public Complex getComplex() {
        return complex;
    }
    public void setLogOdd(){
        ExtractLogOdds lo = new ExtractLogOdds();
        lo.calculateLogOdds();
        double[] logodd = lo.getLo();
        
        //for(int i=0;i<logodd.length;i++)
        //    System.out.print(logodd[i]+";");
        lo.calculateLOTunggal(complex, chain.charAt(0));
    }
    public String setAAI(Aminoacid aa){
        String aaiOfAA = null;
        IO io = new IO();
        AAIndex aai;
        aai= io.extractAAIndexFromFile();
        aai.composeAAI();
        double[] withoutNan =aai.getAaiWithoutNanData(aa.letterToNum());
        aaiOfAA  = aai.printArray1D(withoutNan);
        
        return aaiOfAA;
    }
    public ArrayList<Double>getAAIList(Aminoacid aa){
        IO io = new IO();
        AAIndex aai;
        aai= io.extractAAIndexFromFile();
        aai.composeAAI();
        double[] withoutNan =aai.getAaiWithoutNanData(aa.letterToNum());
        ArrayList<Double> aaid = new ArrayList();
        for(int i=0;i<withoutNan.length;i++){
            aaid.add(withoutNan[i]);
        }
        return aaid;
    }
    public void setSphereExposure() throws IOException{
        ExtractSphereExposure ese;
        ese = new ExtractSphereExposure();
        ArrayList<Aminoacid> listAA = complex.getAminoacidVector();
        for(Aminoacid aa: listAA){
            ese.calculateSE(aa, complex);
        }
        
    }
    public  void setExposedRFromPSAIA(String pdbid, String path, String akhiran){
        ExtractProtusionIndex ePI = new ExtractProtusionIndex();
        ArrayList<PSAI> psailist =ePI.readPIFromFile(this.pdbid,path,akhiran);
        System.out.println("psailist: "+psailist.size());
        ArrayList<PSAI> list =ePI.copyPIOfChain(this.chain, psailist);
        System.out.println("list: "+list.size());
        //untuk setiap aa pada c set rsabytien
        for(int i=0;i<list.size();i++){
            double[]psaiparam =list.get(i).getPsaiParam();
               String residu = this.complex.getAminoacidVector().get(i).getresidueName();
           //    System.out.println("psaiparam"+ psaiparam[0]);
            this.complex.getAminoacidVector().get(i).setRsaTien2013(calculateRSAByTienMSA(residu,psaiparam[0]));
            this.complex.getAminoacidVector().get(i).setRsaRostSander2005(calculateRSARostAndSander(residu,psaiparam[0]));
            this.complex.getAminoacidVector().get(i).setRsaMiller1987(calculateRSAMiller(residu,psaiparam[0]));
            this.complex.getAminoacidVector().get(i).setRsaRose1985(calculateRSARose(residu,psaiparam[0]));
            this.complex.getAminoacidVector().get(i).setsolventAccessibility(psaiparam[0]);
            //set PSAI
            this.complex.getAminoacidVector().get(i).setPsai(list.get(i));
            //if(i<10) {System.out.println("rsa tien:"+ this.complex.getAminoacidVector().get(i).getRsaTien2013());}
        }
                
    }
    public  void setExposedResidueFromPSAIA(){
        ExtractProtusionIndex ePI = new ExtractProtusionIndex();
        ArrayList<PSAI> psailist =ePI.readPIFromFile(this.pdbid);
       // System.out.println("psailist: "+psailist.size());
        ArrayList<PSAI> list =ePI.copyPIOfChain(this.chain, psailist);
        
        //untuk setiap aa pada c set rsabytien
        for(int i=0;i<list.size();i++){
            double[]psaiparam =list.get(i).getPsaiParam();
               String residu = this.complex.getAminoacidVector().get(i).getresidueName();
           //    System.out.println("psaiparam"+ psaiparam[0]);
            this.complex.getAminoacidVector().get(i).setRsaTien2013(calculateRSAByTienMSA(residu,psaiparam[0]));
            this.complex.getAminoacidVector().get(i).setRsaRostSander2005(calculateRSARostAndSander(residu,psaiparam[0]));
            this.complex.getAminoacidVector().get(i).setRsaMiller1987(calculateRSAMiller(residu,psaiparam[0]));
            this.complex.getAminoacidVector().get(i).setRsaRose1985(calculateRSARose(residu,psaiparam[0]));
            this.complex.getAminoacidVector().get(i).setsolventAccessibility(psaiparam[0]);
            //set PSAI
            this.complex.getAminoacidVector().get(i).setPsai(list.get(i));
            //if(i<10) {System.out.println("rsa tien:"+ this.complex.getAminoacidVector().get(i).getRsaTien2013());}
        }
                
    }
    private double calculateRSAMiller(String residueName, double area)
  {
    if ((residueName.equals("ALA")) || (residueName.equals("A"))) {
      return area / 113.0D;
    }
    if ((residueName.equals("ARG")) || (residueName.equals("R"))) {
      return area / 241.0D;
    }
    if ((residueName.equals("ASN")) || (residueName.equals("N"))) {
      return area / 151.0D;
    }
    if ((residueName.equals("ASP")) || (residueName.equals("D"))) {
      return area / 151.0D;
    }
    if ((residueName.equals("CYS")) || (residueName.equals("C"))) {
      return area / 140.0D;
    }
    if ((residueName.equals("GLN")) || (residueName.equals("Q"))) {
      return area / 183.0D;
    }
    if ((residueName.equals("GLU")) || (residueName.equals("E"))) {
      return area / 189.0D;
    }
    if ((residueName.equals("GLY")) || (residueName.equals("G"))) {
      return area / 85.0D;
    }
    if ((residueName.equals("HIS")) || (residueName.equals("H"))) {
      return area / 194.0D;
    }
    if ((residueName.equals("ILE")) || (residueName.equals("I"))) {
      return area / 182.0D;
    }
    if ((residueName.equals("LEU")) || (residueName.equals("L"))) {
      return area / 180.0D;
    }
    if ((residueName.equals("LYS")) || (residueName.equals("K"))) {
      return area / 211.0D;
    }
    if ((residueName.equals("MET")) || (residueName.equals("M"))) {
      return area / 204.0D;
    }
    if ((residueName.equals("PHE")) || (residueName.equals("F"))) {
      return area / 218.0D;
    }
    if ((residueName.equals("PRO")) || (residueName.equals("P"))) {
      return area / 143.0D;
    }
    if ((residueName.equals("SER")) || (residueName.equals("S"))) {
      return area / 122.0D;
    }
    if ((residueName.equals("THR")) || (residueName.equals("T"))) {
      return area / 146.0D;
    }
    if ((residueName.equals("TRP")) || (residueName.equals("W"))) {
      return area / 259.0D;
    }
    if ((residueName.equals("TYR")) || (residueName.equals("Y"))) {
      return area / 229.0D;
    }
    if ((residueName.equals("VAL")) || (residueName.equals("V"))) {
      return area / 160.0D;
    }
    return -1.0D;
  }
  private double calculateRSARose(String residueName, double area)
  {
    if ((residueName.equals("ALA")) || (residueName.equals("A"))) {
      return area / 118.1D;
    }
    if ((residueName.equals("ARG")) || (residueName.equals("R"))) {
      return area / 256.0D;
    }
    if ((residueName.equals("ASN")) || (residueName.equals("N"))) {
      return area / 158.7D;
    }
    if ((residueName.equals("ASP")) || (residueName.equals("D"))) {
      return area / 158.0D;
    }
    if ((residueName.equals("CYS")) || (residueName.equals("C"))) {
      return area / 146.1D;
    }
    if ((residueName.equals("GLN")) || (residueName.equals("Q"))) {
      return area / 186.2D;
    }
    if ((residueName.equals("GLU")) || (residueName.equals("E"))) {
      return area / 193.0D;
    }
    if ((residueName.equals("GLY")) || (residueName.equals("G"))) {
      return area / 88.1D;
    }
    if ((residueName.equals("HIS")) || (residueName.equals("H"))) {
      return area / 202.5D;
    }
    if ((residueName.equals("ILE")) || (residueName.equals("I"))) {
      return area / 181.0D;
    }
    if ((residueName.equals("LEU")) || (residueName.equals("L"))) {
      return area / 193.0D;
    }
    if ((residueName.equals("LYS")) || (residueName.equals("K"))) {
      return area / 225.8D;
    }
    if ((residueName.equals("MET")) || (residueName.equals("M"))) {
      return area / 203.4D;
    }
    if ((residueName.equals("PHE")) || (residueName.equals("F"))) {
      return area / 222.8D;
    }
    if ((residueName.equals("PRO")) || (residueName.equals("P"))) {
      return area / 146.8D;
    }
    if ((residueName.equals("SER")) || (residueName.equals("S"))) {
      return area / 129.8D;
    }
    if ((residueName.equals("THR")) || (residueName.equals("T"))) {
      return area / 152.5D;
    }
    if ((residueName.equals("TRP")) || (residueName.equals("W"))) {
      return area / 266.0D;
    }
    if ((residueName.equals("TYR")) || (residueName.equals("Y"))) {
      return area / 236.0D;
    }
    if ((residueName.equals("VAL")) || (residueName.equals("V"))) {
      return area / 164.5D;
    }
    return -1.0D;
  }
    private double calculateRSAByTienMSA(String residueName, double area)
  {
    if ((residueName.equals("ALA")) || (residueName.equals("A"))) {
      return area / 129.0D;
    }
    if ((residueName.equals("ARG")) || (residueName.equals("R"))) {
      return area / 274.0D;
    }
    if ((residueName.equals("ASN")) || (residueName.equals("N"))) {
      return area / 195.0D;
    }
    if ((residueName.equals("ASP")) || (residueName.equals("D"))) {
      return area / 193.0D;
    }
    if ((residueName.equals("CYS")) || (residueName.equals("C"))) {
      return area / 167.0D;
    }
    if ((residueName.equals("GLN")) || (residueName.equals("Q"))) {
      return area / 223.0D;
    }
    if ((residueName.equals("GLU")) || (residueName.equals("E"))) {
      return area / 225.0D;
    }
    if ((residueName.equals("GLY")) || (residueName.equals("G"))) {
      return area / 104.0D;
    }
    if ((residueName.equals("HIS")) || (residueName.equals("H"))) {
      return area / 224.0D;
    }
    if ((residueName.equals("ILE")) || (residueName.equals("I"))) {
      return area / 197.0D;
    }
    if ((residueName.equals("LEU")) || (residueName.equals("L"))) {
      return area / 201.0D;
    }
    if ((residueName.equals("LYS")) || (residueName.equals("K"))) {
      return area / 236.0D;
    }
    if ((residueName.equals("MET")) || (residueName.equals("M"))) {
      return area / 224.0D;
    }
    if ((residueName.equals("PHE")) || (residueName.equals("F"))) {
      return area / 240.0D;
    }
    if ((residueName.equals("PRO")) || (residueName.equals("P"))) {
      return area / 159.0D;
    }
    if ((residueName.equals("SER")) || (residueName.equals("S"))) {
      return area / 155.0D;
    }
    if ((residueName.equals("THR")) || (residueName.equals("T"))) {
      return area / 172.0D;
    }
    if ((residueName.equals("TRP")) || (residueName.equals("W"))) {
      return area / 285.0D;
    }
    if ((residueName.equals("TYR")) || (residueName.equals("Y"))) {
      return area / 263.0D;
    }
    if ((residueName.equals("VAL")) || (residueName.equals("V"))) {
      return area / 174.0D;
    }
    return -1.0D;
  }
    private double calculateRSARostAndSander(String residueName, double area)
  {
    if ((residueName.equals("ALA")) || (residueName.equals("A"))) {
      return area / 106.0D;
    }
    if ((residueName.equals("ARG")) || (residueName.equals("R"))) {
      return area / 248.0D;
    }
    if ((residueName.equals("ASN")) || (residueName.equals("N"))) {
      return area / 157.0D;
    }
    if ((residueName.equals("ASP")) || (residueName.equals("D"))) {
      return area / 163.0D;
    }
    if ((residueName.equals("CYS")) || (residueName.equals("C"))) {
      return area / 135.0D;
    }
    if ((residueName.equals("GLN")) || (residueName.equals("Q"))) {
      return area / 198.0D;
    }
    if ((residueName.equals("GLU")) || (residueName.equals("E"))) {
      return area / 194.0D;
    }
    if ((residueName.equals("GLY")) || (residueName.equals("G"))) {
      return area / 84.0D;
    }
    if ((residueName.equals("HIS")) || (residueName.equals("H"))) {
      return area / 184.0D;
    }
    if ((residueName.equals("ILE")) || (residueName.equals("I"))) {
      return area / 169.0D;
    }
    if ((residueName.equals("LEU")) || (residueName.equals("L"))) {
      return area / 164.0D;
    }
    if ((residueName.equals("LYS")) || (residueName.equals("K"))) {
      return area / 205.0D;
    }
    if ((residueName.equals("MET")) || (residueName.equals("M"))) {
      return area / 188.0D;
    }
    if ((residueName.equals("PHE")) || (residueName.equals("F"))) {
      return area / 197.0D;
    }
    if ((residueName.equals("PRO")) || (residueName.equals("P"))) {
      return area / 136.0D;
    }
    if ((residueName.equals("SER")) || (residueName.equals("S"))) {
      return area / 130.0D;
    }
    if ((residueName.equals("THR")) || (residueName.equals("T"))) {
      return area / 142.0D;
    }
    if ((residueName.equals("TRP")) || (residueName.equals("W"))) {
      return area / 227.0D;
    }
    if ((residueName.equals("TYR")) || (residueName.equals("Y"))) {
      return area / 222.0D;
    }
    if ((residueName.equals("VAL")) || (residueName.equals("V"))) {
      return area / 142.0D;
    }
    return -1.0D;
  }
}
