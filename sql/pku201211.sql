SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='TRADITIONAL,ALLOW_INVALID_DATES';

CREATE SCHEMA IF NOT EXISTS `pkudb` DEFAULT CHARACTER SET utf8 COLLATE utf8_general_ci ;
SHOW WARNINGS;
USE `pkudb` ;

-- -----------------------------------------------------
-- Table `pkudb`.`Samples`
-- -----------------------------------------------------
CREATE  TABLE IF NOT EXISTS `pkudb`.`Samples` (
  `SampleID` CHAR(12) NOT NULL ,
  `SourceID` CHAR(4) NULL ,
  `SourceCnt` SMALLINT(4) NULL ,
  `TissueID` CHAR(2) NULL ,
  `TissueCnt` TINYINT(2) NULL ,
  `DateCollected` DATETIME NULL ,
  `DateReceived` DATETIME NULL ,
  `DateFrozen` DATETIME NULL ,
  `SampleCodePre` VARCHAR(45) NULL ,
  `TissueReal` VARCHAR(45) NULL ,
  `AnimalID` VARCHAR(45) NULL COMMENT 'Def = SourceID . SourceCnt' ,
  `OldID` VARCHAR(45) NULL ,
  `LabelPrintedCnt` INT NOT NULL DEFAULT 0 ,
  `ivFreezer` VARCHAR(45) NULL ,
  `ivShelf` VARCHAR(45) NULL ,
  `ivRack` VARCHAR(45) NULL ,
  `ivBox` VARCHAR(45) NULL ,
  `ivPosition` VARCHAR(45) NULL ,
  `Description` VARCHAR(4095) NULL ,
  `ProjectName` VARCHAR(45) NULL ,
  `Draw` VARCHAR(45) NULL ,
  `Vial` VARCHAR(45) NULL ,
  `TransportTemperature` DECIMAL(4,1) NULL ,
  `Quality` VARCHAR(45) NULL ,
  `VolumeInTubeUL` DECIMAL(5,1) NULL ,
  `CreatedBy` VARCHAR(45) NULL ,
  `CreatedAt` DATETIME NULL ,
  `ModifiedBy` VARCHAR(45) NULL ,
  `ModifiedAt` DATETIME NULL ,
  `Depleted` TINYINT(1) NULL ,
  PRIMARY KEY (`SampleID`) ,
  UNIQUE INDEX `SampleCode_UNIQUE` (`SampleID` ASC) )
ENGINE = InnoDB;

SHOW WARNINGS;

-- -----------------------------------------------------
-- Table `pkudb`.`AdditionalSampleInfo`
-- -----------------------------------------------------
CREATE  TABLE IF NOT EXISTS `pkudb`.`AdditionalSampleInfo` (
  `UID` INT NOT NULL AUTO_INCREMENT ,
  `SampleID` CHAR(12) NOT NULL ,
  `Key` VARCHAR(45) NOT NULL ,
  `Value` VARCHAR(4095) NULL ,
  UNIQUE INDEX `SampleID_UNIQUE` (`SampleID` ASC) ,
  PRIMARY KEY (`UID`) ,
  UNIQUE INDEX `UID_UNIQUE` (`UID` ASC) )
ENGINE = InnoDB;

SHOW WARNINGS;

-- -----------------------------------------------------
-- Table `pkudb`.`Project`
-- -----------------------------------------------------
CREATE  TABLE IF NOT EXISTS `pkudb`.`Project` (
  `UID` INT NOT NULL ,
  `ProjectName` VARCHAR(45) NOT NULL ,
  `Started` DATETIME NULL ,
  `Finished` DATETIME NULL ,
  `Status` VARCHAR(45) NULL ,
  `Comments` VARCHAR(4095) NULL ,
  PRIMARY KEY (`UID`) ,
  UNIQUE INDEX `UID_UNIQUE` (`UID` ASC) ,
  UNIQUE INDEX `ProjectName_UNIQUE` (`ProjectName` ASC) )
ENGINE = InnoDB;

SHOW WARNINGS;

-- -----------------------------------------------------
-- Table `pkudb`.`SampleRemoval`
-- -----------------------------------------------------
CREATE  TABLE IF NOT EXISTS `pkudb`.`SampleRemoval` (
  `SampleID` CHAR(12) NOT NULL ,
  `rmBy` VARCHAR(45) NULL ,
  `rmAt` VARCHAR(45) NULL ,
  `rmReason` VARCHAR(45) NULL ,
  `rmAmount` VARCHAR(45) NULL ,
  `Comments` VARCHAR(4095) NULL ,
  PRIMARY KEY (`SampleID`) ,
  UNIQUE INDEX `SampleID_UNIQUE` (`SampleID` ASC) )
ENGINE = InnoDB;

SHOW WARNINGS;

-- -----------------------------------------------------
-- Table `pkudb`.`Contacts`
-- -----------------------------------------------------
CREATE  TABLE IF NOT EXISTS `pkudb`.`Contacts` (
  `SampleID` CHAR(12) NOT NULL ,
  `LabPerson` VARCHAR(45) NOT NULL ,
  `Collector` VARCHAR(45) NULL ,
  `Information` VARCHAR(4095) NULL ,
  PRIMARY KEY (`SampleID`) ,
  UNIQUE INDEX `SampleID_UNIQUE` (`SampleID` ASC) )
ENGINE = InnoDB;

SHOW WARNINGS;

-- -----------------------------------------------------
-- Table `pkudb`.`Geography`
-- -----------------------------------------------------
CREATE  TABLE IF NOT EXISTS `pkudb`.`Geography` (
  `SampleID` CHAR(12) NOT NULL ,
  `Continent` VARCHAR(45) NULL ,
  `Country` VARCHAR(45) NULL ,
  `StateProvince` VARCHAR(45) NULL ,
  `Region` VARCHAR(45) NULL ,
  `SubRegion` VARCHAR(45) NULL ,
  `Info` VARCHAR(45) NULL ,
  `Comments` VARCHAR(4095) NULL ,
  `Latitude` VARCHAR(45) NULL ,
  `Longitude` VARCHAR(45) NULL ,
  `Elevation` VARCHAR(45) NULL ,
  PRIMARY KEY (`SampleID`) ,
  UNIQUE INDEX `SampleID_UNIQUE` (`SampleID` ASC) )
ENGINE = InnoDB;

SHOW WARNINGS;

-- -----------------------------------------------------
-- Table `pkudb`.`ProjectList`
-- -----------------------------------------------------
CREATE  TABLE IF NOT EXISTS `pkudb`.`ProjectList` (
  `ProjectID` INT NOT NULL ,
  `SampleID` CHAR(12) NULL ,
  `AnimalID` VARCHAR(45) NULL ,
  INDEX `UsedSamples` USING HASH (`SampleID` ASC, `AnimalID` ASC) )
ENGINE = InnoDB;

SHOW WARNINGS;

-- -----------------------------------------------------
-- Table `pkudb`.`TissueCode`
-- -----------------------------------------------------
CREATE  TABLE IF NOT EXISTS `pkudb`.`TissueCode` (
  `TissueID` CHAR(2) NOT NULL ,
  `Tissue` VARCHAR(45) NULL ,
  `Description` VARCHAR(4095) NULL ,
  PRIMARY KEY (`TissueID`) ,
  UNIQUE INDEX `TissueID_UNIQUE` (`TissueID` ASC) ,
  UNIQUE INDEX `Tissue_UNIQUE` (`Tissue` ASC) )
ENGINE = InnoDB;

SHOW WARNINGS;

-- -----------------------------------------------------
-- Table `pkudb`.`Animals`
-- -----------------------------------------------------
CREATE  TABLE IF NOT EXISTS `pkudb`.`Animals` (
  `AnimalID` VARCHAR(45) NOT NULL ,
  `Sex` CHAR(1) NULL ,
  `SexDeterminedBy` VARCHAR(45) NULL ,
  `SexComments` VARCHAR(4095) NULL ,
  `DateBirth` DATETIME NULL ,
  `DateDeath` DATETIME NULL ,
  `Date1stSampling` DATETIME NULL ,
  `CauseOfDeath` VARCHAR(4095) NULL ,
  `Comments` VARCHAR(4095) NULL ,
  `StatusBirth` VARCHAR(45) NULL ,
  `StatusCurrent` VARCHAR(45) NULL ,
  PRIMARY KEY (`AnimalID`) ,
  UNIQUE INDEX `AnimalID_UNIQUE` (`AnimalID` ASC) )
ENGINE = InnoDB;

SHOW WARNINGS;

-- -----------------------------------------------------
-- Table `pkudb`.`OtherID`
-- -----------------------------------------------------
CREATE  TABLE IF NOT EXISTS `pkudb`.`OtherID` (
  `OID` VARCHAR(45) NOT NULL ,
  `identifier` VARCHAR(45) NULL ,
  `Comments` VARCHAR(45) NULL ,
  UNIQUE INDEX `OID_UNIQUE` (`OID` ASC) )
ENGINE = InnoDB;

SHOW WARNINGS;

-- -----------------------------------------------------
-- Table `pkudb`.`OtherIDList`
-- -----------------------------------------------------
CREATE  TABLE IF NOT EXISTS `pkudb`.`OtherIDList` (
  `OID` VARCHAR(45) NOT NULL ,
  `SampleID` CHAR(12) NULL ,
  `AnimalID` VARCHAR(45) NULL ,
  INDEX `UsedSamples` USING HASH (`SampleID` ASC, `AnimalID` ASC) ,
  UNIQUE INDEX `OID_UNIQUE` (`OID` ASC) ,
  INDEX `toSamples_idx` (`SampleID` ASC) ,
  INDEX `toAnimals_idx` (`AnimalID` ASC) ,
  CONSTRAINT `toSamples`
    FOREIGN KEY (`SampleID` )
    REFERENCES `pkudb`.`Samples` (`SampleID` )
    ON DELETE NO ACTION
    ON UPDATE CASCADE,
  CONSTRAINT `toAnimals`
    FOREIGN KEY (`AnimalID` )
    REFERENCES `pkudb`.`Animals` (`AnimalID` )
    ON DELETE NO ACTION
    ON UPDATE CASCADE)
ENGINE = InnoDB;

SHOW WARNINGS;
USE `pkudb`;

DELIMITER $$
SHOW WARNINGS$$
USE `pkudb`$$


CREATE TRIGGER SampleID_Split BEFORE INSERT ON Samples
FOR EACH ROW
BEGIN
	SET NEW.SampleID=UPPER(NEW.SampleID);
	SET NEW.CreatedAt = CURRENT_TIMESTAMP;
	IF NEW.SourceID IS NULL OR NEW.SourceID='' THEN
		SET NEW.SourceID=SUBSTRING(NEW.SampleID,1,4),
			NEW.SourceCnt=SUBSTRING(NEW.SampleID,5,4),
			NEW.TissueID=SUBSTRING(NEW.SampleID,9,2),
			NEW.TissueCnt=SUBSTRING(NEW.SampleID,11);
	END IF;
	IF NEW.AnimalID IS NULL OR NEW.AnimalID='' THEN
		SET NEW.AnimalID=CONCAT(NEW.SourceID,LPAD(NEW.SourceCnt,4,'0'));
	END IF;
	INSERT IGNORE INTO Animals SET AnimalID=NEW.AnimalID;
	INSERT IGNORE INTO TissueCode SET TissueID=NEW.TissueID;
END; $$

SHOW WARNINGS$$

CREATE TRIGGER Animals_Date1stSampling AFTER INSERT ON Samples
FOR EACH ROW
BEGIN
	DECLARE Oldata DATETIME;
	IF NEW.DateCollected IS NOT NULL THEN
		SELECT Date1stSampling FROM Animals WHERE AnimalID=NEW.AnimalID INTO Oldata;
		IF Oldata > NEW.DateCollected THEN
			UPDATE IGNORE Animals SET Date1stSampling=NEW.DateCollected WHERE AnimalID=NEW.AnimalID;
		END IF;
	END IF;
END; $$

DELIMITER ;


SET SQL_MODE=@OLD_SQL_MODE;
SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;

-- -----------------------------------------------------
-- Test Data
-- -----------------------------------------------------
USE `pkudb`;
INSERT INTO Samples (SampleID) VALUES ('gdxj0000bl00');
INSERT INTO Samples (SampleID) VALUES ('gdxj9999ms99');
INSERT INTO Samples (SampleID) VALUES ('gdxj0001bl01');
INSERT INTO Samples (SampleID,DateCollected) VALUES ('gdxj0001bl02','2012-07-05 18:11:12');
SELECT AnimalID,Date1stSampling FROM Animals;
INSERT INTO Samples (SampleID,DateCollected) VALUES ('gdxj0001bl12','2012-06-05');
SELECT AnimalID,Date1stSampling FROM Animals;
INSERT INTO Samples (SampleID,DateCollected) VALUES ('gdxj0001bl32','2012-08-05 17:11:12');
SELECT AnimalID,Date1stSampling FROM Animals;
SELECT SampleID,SourceID,SourceCnt,TissueID,TissueCnt,AnimalID,DateCollected,LabelPrinted FROM Samples;
