package ComplexStructure;
public class Atom
{
  private String atomName;
  private String residueName;
  private char chainName;
  private String residueID;
  private double bFactor;
  private double x;
  private double y;
  private double z;
  private int type;
  
  public void setatomName(String atomName)
  {
    this.atomName = atomName;
  }
  
  public void setbFactor(double bFactor)
  {
    this.bFactor = bFactor;
  }
  
  public void setresidueName(String residueName)
  {
    this.residueName = residueName;
  }
  
  public void setchainName(char chainName)
  {
    this.chainName = chainName;
  }
  
  public void setresidueID(String residueID)
  {
    this.residueID = residueID;
  }
  
  public int getType()
  {
    return this.type;
  }
  
  public void setType(int type)
  {
    this.type = type;
  }
  
  public String getatomName()
  {
    return this.atomName;
  }
  
  public String getresidueName()
  {
    return this.residueName;
  }
  
  public char getchainID()
  {
    return this.chainName;
  }
  
  public String getresidueID()
  {
    return this.residueID;
  }
  
  public double getbFactor()
  {
    return this.bFactor;
  }
  
  public double getX()
  {
    return this.x;
  }
  
  public void setX(double x)
  {
    this.x = x;
  }
  
  public double getY()
  {
    return this.y;
  }
  
  public void setY(double y)
  {
    this.y = y;
  }
  
  public double getZ()
  {
    return this.z;
  }
  
  public void setZ(double z)
  {
    this.z = z;
  }
}